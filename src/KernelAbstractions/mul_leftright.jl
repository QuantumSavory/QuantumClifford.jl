
#=============================================================================#
# CAUTION: Keep in mind that mutable = order_right_left ? right : left.
# TODO: Make the parameters keyword arguments once support becomes available.
KA.@kernel inbounds = true unsafe_indices = true function kernel_mul!(
	mutable_phases, mutable_xzs, @Const(const_xzs), @Const(stride),
	::Val{order_right_left}, ::Val{phases},
	::Val{block_size}, ::Val{batch_size}
	) where {order_right_left, phases, block_size, batch_size}

	# unsafe_indices is required for shared_memory_reduce, do this manually.
	global_position = global_index(
		KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
		)
	start = global_position[1]
	j_mutable = global_position[2]
	j_const = 1
	count = KA.@uniform (size(const_xzs, 1) >> 1)
	if phases
		low = KA.@uniform (zero(eltype(const_xzs)))
		high = KA.@uniform (zero(eltype(const_xzs)))
	end

	for (i, _) in zip(start : stride : count, one(batch_size) : batch_size)
		if order_right_left
			x_left = const_xzs[i, j_const]
			z_left = const_xzs[i + count, j_const]
			x_right = mutable_xzs[i, j_mutable]
			z_right = mutable_xzs[i + count, j_mutable]
		else
			x_left = mutable_xzs[i, j_mutable]
			z_left = mutable_xzs[i + count, j_mutable]
			x_right = const_xzs[i, j_const]
			z_right = const_xzs[i + count, j_const]
		end

		x_new = xor(x_left, x_right)
		z_new = xor(z_left, z_right)
		if phases
			xl_zr = x_left & z_right
			merged = xor(xl_zr, z_left & x_right)
			high =
				xor(high, xor(low, x_new, z_new, xl_zr) & merged)
			low = xor(low, merged)
		end

		mutable_xzs[i, j_mutable] = x_new
		mutable_xzs[i + count, j_mutable] = z_new
	end

	if phases
		value::DeviceUnsigned = (count_ones(high) << 1) + count_ones(low)
		buffer = KA.@localmem DeviceUnsigned (block_size,)
		index = KA.@index(Local, Linear)
		value = shared_memory_reduce!(+, buffer, value, index, Val(block_size))
		if index == one(index)
			Atomix.@atomic mutable_phases[j_mutable] += value
		end
	end

end

function device_mul!(
	mutable_phases::AbstractGPUArray{<: Unsigned},
	mutable_xzs::AbstractGPUArray{T},
	const_phases::AbstractGPUArray{<: Unsigned},
	const_xzs::AbstractGPUArray{T};
	order_right_left::Val{right_left}, phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	)::Nothing where {T <: Unsigned, right_left, phase_B, block_SZ, batch_SZ}

	backend = KA.get_backend(mutable_xzs)
	dim_x = max(block_SZ, cld(size(mutable_xzs, 1) >> 1, batch_SZ))
	dim_y = size(mutable_xzs, 2)
	# This resolves a race condition if the phases alias each other.
	if phase_B
		const_phases = copy(const_phases)
	end
	mul! = kernel_mul!(backend, (block_SZ, 1))
	mul!(
		mutable_phases, mutable_xzs, const_xzs, dim_x,
		order_right_left, phases, block_size, batch_size;
		workgroupsize = (block_SZ, 1), ndrange = (dim_x, dim_y)
		)
	if phase_B
		transform! = kernel_transform!(backend)
		transform!(mod_4_sum!, mutable_phases, const_phases; ndrange = dim_y)
	end
	return nothing

end

# CAUTION: Meta-programming is utilised to combine left/right definitions.
for (direction, right_left) in ((:left!, true), (:right!, false))

# Avoid creating messy variables in the outermost scope.
local sym = :mul_
local prefix_sym = :do_mul_

#==============================================================================
RETURNS PAULI OPERATOR
==============================================================================#

# PauliOperator - PauliOperator
@eval @inline function $(Symbol(sym, direction))(
	u::DevicePauliOperator, v::DevicePauliOperator;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	u.nqubits == v.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	device_mul!(
		u.phase, u.xz, v.phase, v.xz;
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DevicePauliOperator, v::DevicePauliOperator;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	device_mul!(
		u.phase, u.xz, v.phase, v.xz;
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

# PauliOperator - Tableau[i]
@eval @inline function $(Symbol(sym, direction))(
	u::DevicePauliOperator, v::DeviceTableau, i;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	u.nqubits == v.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	device_mul!(
		u.phase, u.xz, (@view v.phases[i]), (@view v.xzs[:, i]);
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DevicePauliOperator, v::DeviceTableau, i;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	@inbounds device_mul!(
		u.phase, u.xz, (@view v.phases[i]), (@view v.xzs[:, i]);
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

# PauliOperator - AbstractStabilizer[i]
@eval @inline function $(Symbol(sym, direction))(
	u::DevicePauliOperator, v::DeviceAbstractStabilizer, i;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	$(Symbol(sym, direction))(
		u, v.tab, i;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DevicePauliOperator, v::DeviceAbstractStabilizer, i;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	return $(Symbol(prefix_sym, direction))(
		u, v.tab, i;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

#==============================================================================
RETURNS TABLEAU
==============================================================================#

# Tableau - PauliOperator
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceTableau, v::DevicePauliOperator;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	u.nqubits == v.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	device_mul!(
		u.phases, u.xzs, v.phase, v.xz;
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceTableau, v::DevicePauliOperator;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	device_mul!(
		u.phases, u.xzs, v.phase, v.xz;
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

# Tableau[m] - Tableau[n]
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceTableau, m, v::DeviceTableau, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	u.nqubits == v.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	device_mul!(
		(@view u.phases[m]), (@view u.xzs[:, m]),
		(@view v.phases[n]), (@view v.xzs[:, n]);
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceTableau, m, v::DeviceTableau, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	@inbounds device_mul!(
		(@view u.phases[m]), (@view u.xzs[:, m]),
		(@view v.phases[n]), (@view v.xzs[:, n]);
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

# Tableau[m] - AbstractStabilizer[n]
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceTableau, m, v::DeviceAbstractStabilizer, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	return $(Symbol(sym, direction))(
		u, m, v.tab, n;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceTableau, m, v::DeviceAbstractStabilizer, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	return $(Symbol(prefix_sym, direction))(
		u, m, v.tab, n;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

# Tableau[m] - Self[n]
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceTableau, m, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}


	device_mul!(
		(@view u.phases[m]), (@view u.xzs[:, m]),
		(@view u.phases[n]), (@view u.xzs[:, n]);
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceTableau, m, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	@inbounds device_mul!(
		(@view u.phases[m]), (@view u.xzs[:, m]),
		(@view u.phases[n]), (@view u.xzs[:, n]);
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

#==============================================================================
RETURNS (MIXED) (DE)STABILIZER
==============================================================================#

# AbstractStabilizer - PauliOperator
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceAbstractStabilizer, v::DevicePauliOperator;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}


	$(Symbol(sym, direction))(
		u.tab, v;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceAbstractStabilizer, v::DevicePauliOperator;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}


	$(Symbol(prefix_sym, direction))(
		u.tab, v;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	return u

end

# AbstractStabilizer[m] - Tableau[n]
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceAbstractStabilizer, m, v::DeviceTableau, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}


	$(Symbol(sym, direction))(
		u.tab, m, v, n;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceAbstractStabilizer, m, v::DeviceTableau, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}


	$(Symbol(prefix_sym, direction))(
		u.tab, m, v, n;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	return u

end

# AbstractStabilizer[m] - AbstractStabilizer[n]
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceAbstractStabilizer, m, v::DeviceAbstractStabilizer, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	$(Symbol(sym, direction))(
		u.tab, m, v.tab, n;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceAbstractStabilizer, m, v::DeviceAbstractStabilizer, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	$(Symbol(prefix_sym, direction))(
		u.tab, m, v.tab, n;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	return u

end

#==============================================================================
RETURNS (MIXED) STABILIZER
==============================================================================#

# (Mixed)Stabilizer[m] - Self[n]
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceUnionStabilizer, m, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}


	$(Symbol(sym, direction))(
		u.tab, m, n;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceUnionStabilizer, m, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	$(Symbol(prefix_sym, direction))(
		u.tab, m, n;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	return u

end

#==============================================================================
RETURNS (MIXED) DESTABILIZER
==============================================================================#

# (Mixed)Destabilizer[m] - Self[n]
@eval @inline function $(Symbol(sym, direction))(
	u::DeviceUnionDestabilizer, m, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	p, nqubits, xzs = u.tab.phases, u.tab.nqubits, u.tab.xzs
	# Swapping the order of the indices is intentional.
	a_phase, a_xz = (@view p[n]), (@view xzs[:, n])
	b_phase, b_xz = (@view p[m]), (@view xzs[:, m])
	# This trick enables supporting CartesianIndex.
	p, xzs = (@view p[nqubits + 1 : end]), (@view xzs[nqubits + 1 : end])
	c_phase, c_xz = (@view p[m]), (@view xzs[:, m])
	d_phase, d_xz = (@view p[n]), (@view xzs[:, n])
	device_mul!(
		a_phase, a_xz, b_phase, b_xz;
		order_right_left = Val($right_left), phases = Val(false),
		block_size = block_size, batch_size = batch_size
		)
	device_mul!(
		c_phase, c_xz, d_phase, d_xz;
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

@eval @inline function $(Symbol(prefix_sym, direction))(
	u::DeviceUnionDestabilizer, m, n;
	phases::Val{phase_B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {phase_B, block_SZ, batch_SZ}

	p, nqubits, xzs = u.tab.phases, u.tab.nqubits, u.tab.xzs
	# Swapping the order of the indices is intentional.
	@inbounds device_mul!(
		(@view p[n]), (@view xzs[:, n]), (@view p[m]), (@view xzs[:, m]);
		order_right_left = Val($right_left), phases = Val(false),
		block_size = block_size, batch_size = batch_size
		)
	# This trick enables supporting CartesianIndex.
	@inbounds p, xzs =
		(@view p[nqubits + 1 : end]), (@view xzs[nqubits + 1 : end])
	@inbounds device_mul!(
		(@view p[m]), (@view xzs[:, m]), (@view p[n]), (@view xzs[:, n]);
		order_right_left = Val($right_left), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return u

end

# Marks the end for (direction, right_left)
end
#=============================================================================#
