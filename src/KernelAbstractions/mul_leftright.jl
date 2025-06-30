
#=============================================================================#
# CAUTION: Keep in mind that mutable = order_right_left ? right : left.
KA.@kernel inbounds = true unsafe_indices = true function kernel_mul!(
	mutable_phases, mutable_xzs, @Const(const_xzs), @Const(stride),
	::Val{order_right_left}, ::Val{phases},
	::Val{block_size}, ::Val{batch_size}
	) where {order_right_left, block_size, phases, batch_size}

	# unsafe_indices is required for shared_memory_reduce, do this manually.
	global_pos = global_index(
		KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
		)
	start = global_pos[1]
	j_mutable = global_pos[2]
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
		value = shared_memory_reduce!(
			+, buffer, value, index, Val(block_size)
			)
		if index == one(index)
			Atomix.@atomic mutable_phases[j_mutable] +=
				convert(eltype(mutable_phases), value)
		end
	end

end

function mul_device!(
	mutable_phases::AbstractArray{<: Unsigned}, mutable_xzs::AbstractArray{T},
	const_phases::AbstractArray{<: Unsigned}, const_xzs::AbstractArray{T};
	order_right_left::Val{RL}, phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	)::Nothing where {T <: Unsigned, RL, B, block_SZ, batch_SZ}

	backend = KA.get_backend(mutable_xzs)
	dim_x = max(block_SZ, cld(size(mutable_xzs, 1) >> 1, batch_SZ))
	dim_y = size(mutable_xzs, 2)
	# This resolves a race condition if the phases alias each other.
	if B
		const_phases = copy(const_phases)
	end
	mul! = kernel_mul!(backend, (block_SZ, 1))
	mul!(
		mutable_phases, mutable_xzs, const_xzs, dim_x,
		order_right_left, phases, block_size, batch_size;
		workgroupsize = (block_SZ, 1), ndrange = (dim_x, dim_y)
		)
	if B
		transform! = kernel_transform!(backend)
		transform!(mod_4_sum!, mutable_phases, const_phases; ndrange = dim_y)
	end
	return nothing

end

#==============================================================================
							RETURNS PAULI OPERATOR
==============================================================================#

@inline function mul_left!(
	r::PauliOperator{R_P, R_XZ}, l::PauliOperator{L_P, L_XZ};
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		R_P <: AbstractGPUArray, R_XZ <: AbstractGPUArray,
		L_P <: AbstractGPUArray, L_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	l.nqubits == r.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	return do_mul_left!(
		r, l;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_left!(
	r::PauliOperator{R_P, R_XZ}, l::PauliOperator{L_P, L_XZ};
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		R_P <: AbstractGPUArray, R_XZ <: AbstractGPUArray,
		L_P <: AbstractGPUArray, L_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		r.phase, r.xz, l.phase, l.xz;
		order_right_left = Val(true), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return r

end

@inline function mul_right!(
	l::PauliOperator{L_P, L_XZ}, r::PauliOperator{R_P, R_XZ};
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		L_P <: AbstractGPUArray, L_XZ <: AbstractGPUArray,
		R_P <: AbstractGPUArray, R_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	l.nqubits == r.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	return do_mul_right!(
		l, r;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_right!(
	l::PauliOperator{L_P, L_XZ}, r::PauliOperator{R_P, R_XZ};
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		L_P <: AbstractGPUArray, L_XZ <: AbstractGPUArray,
		R_P <: AbstractGPUArray, R_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		l.phase, l.xz, r.phase, r.xz;
		order_right_left = Val(false), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return l

end

@inline function mul_left!(
	r::PauliOperator{R_P, R_XZ}, l::Tableau{L_P, L_XZ}, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		R_P <: AbstractGPUArray, R_XZ <: AbstractGPUArray,
		L_P <: AbstractGPUArray, L_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	l.nqubits == r.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	return do_mul_left!(
		r, l, i;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_left!(
	r::PauliOperator{R_P, R_XZ}, l::Tableau{L_P, L_XZ}, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		R_P <: AbstractGPUArray, R_XZ <: AbstractGPUArray,
		L_P <: AbstractGPUArray, L_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		r.phase, r.xz, (@view l.phases[i]), (@view l.xzs[:, i]);
		order_right_left = Val(true), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return r

end

@inline function mul_right!(
	l::PauliOperator{L_P, L_XZ}, r::Tableau{R_P, R_XZ}, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		L_P <: AbstractGPUArray, L_XZ <: AbstractGPUArray,
		R_P <: AbstractGPUArray, R_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	l.nqubits == r.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	return do_mul_right!(
		l, r, i;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_right!(
	l::PauliOperator{L_P, L_XZ}, r::Tableau{R_P, R_XZ}, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		L_P <: AbstractGPUArray, L_XZ <: AbstractGPUArray,
		R_P <: AbstractGPUArray, R_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		l.phase, l.xz, (@view r.phases[i]), (@view r.xzs[:, i]);
		order_right_left = Val(false), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return l

end

#==============================================================================
								RETURNS TABLEAU
==============================================================================#

@inline function mul_left!(
	s::Tableau{S_P, S_XZ}, m, t::Tableau{T_P, T_XZ}, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	s.nqubits == t.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	return do_mul_left!(
		s, m, t, i;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_left!(
	s::Tableau{S_P, S_XZ}, m, t::Tableau{T_P, T_XZ}, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		(@view s.phases[m]), (@view s.xzs[:, m]),
		(@view t.phases[i]), (@view t.xzs[:, i]);
		order_right_left = Val(true), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return s

end

@inline function mul_right!(
	s::Tableau{S_P, S_XZ}, m, t::Tableau{T_P, T_XZ}, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	s.nqubits == t.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	return do_mul_right!(
		s, m, t, i;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_right!(
	s::Tableau{S_P, S_XZ}, m, t::Tableau{T_P, T_XZ}, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		(@view s.phases[m]), (@view s.xzs[:, m]),
		(@view t.phases[i]), (@view t.xzs[:, i]);
		order_right_left = Val(false), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return s

end

@inline function mul_left!(
	s::Tableau{S_P, S_XZ}, m, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	return do_mul_left!(
		s, m, t, i;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_left!(
	s::Tableau{S_P, S_XZ}, m, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		(@view s.phases[m]), (@view s.xzs[:, m]),
		(@view s.phases[i]), (@view s.xzs[:, i]);
		order_right_left = Val(true), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return s

end

@inline function mul_right!(
	s::Tableau{S_P, S_XZ}, m, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	return do_mul_right!(
		s, m, t, i;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_right!(
	s::Tableau{S_P, S_XZ}, m, i;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		(@view s.phases[m]), (@view s.xzs[:, m]),
		(@view s.phases[i]), (@view s.xzs[:, i]);
		order_right_left = Val(false), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return s

end

@inline function mul_left!(
	s::Tableau{S_P, S_XZ}, p::PauliOperator{P_P, P_XZ};
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		P_P <: AbstractGPUArray, P_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	s.nqubits == p.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	return do_mul_left!(
		s, p;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_left!(
	s::Tableau{S_P, S_XZ}, p::PauliOperator{P_P, P_XZ};
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		P_P <: AbstractGPUArray, P_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		s.phases, s.xzs, p.phase, p.xz;
		order_right_left = Val(true), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return s

end

@inline function mul_right!(
	s::Tableau{S_P, S_XZ}, p::PauliOperator{P_P, P_XZ};
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		P_P <: AbstractGPUArray, P_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	s.nqubits == p.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
	return do_mul_right!(
		s, p;
		phases = phases, block_size = block_size, batch_size = batch_size
		)

end

@inline function do_mul_right!(
	s::Tableau{S_P, S_XZ}, p::PauliOperator{P_P, P_XZ};
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(default_block_size),
	batch_size::Val{batch_SZ} = Val(default_batch_size)
	) where {
		S_P <: AbstractGPUArray, S_XZ <: AbstractGPUArray,
		P_P <: AbstractGPUArray, P_XZ <: AbstractGPUArray,
		B, block_SZ, batch_SZ
		}

	mul_device!(
		s.phases, s.xzs, p.phase, p.xz;
		order_right_left = Val(false), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	return s

end
#=============================================================================#
