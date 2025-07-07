
#=============================================================================#
# For use whenever KA.@index(Global, NTuple) is unavailable.
# At the moment, JET always complains unless unsafe_indices = true is set.
@inline global_index(block_index, block_dim, thread_index) =
	(block_index .- one(eltype(block_index))) .* block_dim .+ thread_index

# This would have been called "map" but that name is already reserved.
KA.@kernel inbounds = true unsafe_indices = true function kernel_transform!(
	f!, target, @Const(auxiliary)
	)

	index = global_index(
		KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
		)
	f!(target, auxiliary, index)

end

# CAUTION: Requires block_size == length(buffer) == prod(KA.@groupsize()).
# CAUTION: Requires unsafe_indices = true if num_active_threads < block_size.
# TODO: Overhaul once __ctx__ is no longer necessary for runtime queries.
# TODO: Revisit once warp level primitives are supported.
@inline function shared_memory_reduce!(
	f, buffer::AbstractArray{T}, value::T, index, ::Val{block_size}
	) where {T, block_size}

	@inbounds buffer[index] = value

	# This branch is a power of 2, take the quick route.
	if count_ones(block_size) == one(block_size)

	# This is messy but only for-loop unrolling is supported.
	KA.Extras.@unroll for bit = one(block_size) : trailing_zeros(block_size)
		stride = KA.@uniform ((2)^(trailing_zeros(block_size) - bit))
		# The call to KA.@synchronize is ALL or NOTHING.
		KA.@synchronize()
		if index <= stride
			@inbounds value = f(value, buffer[index + stride])
			@inbounds buffer[index] = value
		end
	end

	else

	current = KA.@uniform (block_size)
	# TODO: Unroll this branch should it become possible.
	while current > one(current)
		# The call to KA.@synchronize is ALL or NOTHING.
		KA.@synchronize()
		# This splits into even/odd steps. Should target SALU.
		if iseven(current)
			current = KA.@uniform (current >> one(current))
			if index <= current
				@inbounds value = f(value, buffer[index + current])
				@inbounds buffer[index] = value
			end
		else
			current = KA.@uniform ((current >> one(current)) + one(current))
			# The strict inequality is intentional.
			if index < current
				@inbounds value = f(value, buffer[index + current])
				@inbounds buffer[index] = value
			end
		end
	end

	end
	return value

end

#==============================================================================
							COMMON COMPUTATIONS
==============================================================================#

# Anonymous functions trigger recompilation. Hence, Separate them out.
@inline function mod_4_sum!(t, a, pos)
	i = pos[1]
	if i <= length(t)
		j = length(a) > 1 ? i : 1
		@inbounds t[i] = (t[i] + a[j]) & 0x3
	end
end
#=============================================================================#
