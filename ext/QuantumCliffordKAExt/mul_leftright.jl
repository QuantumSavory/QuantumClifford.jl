
#=============================================================================#
# TODO: include the unsafe functions once the main package establishes them.
import QuantumClifford: mul_left!, mul_right!

# CAUTION: Keep in mind that mutable = order_right_left ? right : left
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
    j_mutable = ifelse(size(mutable_xzs, 2) > 1, global_position[2], 1)
    j_const = ifelse(size(const_xzs, 2) > 1, global_position[2], 1)
    count = KA.@uniform (size(mutable_xzs, 1) >> 1)
    if phases
        low = KA.@uniform (zero(eltype(mutable_xzs)))
        high = KA.@uniform (zero(eltype(mutable_xzs)))
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
            @atomic mutable_phases[j_mutable] += value
        end
    end

end

# CAUTION: Requires either rows(const) == 1 or rows(const) == rows(mutable)
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
    if phase_B
        transform! = kernel_transform!(backend)
        transform!(
            mod_4_sum!, mutable_phases, const_phases; ndrange = dim_y
            )
    end
    mul! = kernel_mul!(backend, (block_SZ, 1))
    mul!(
        mutable_phases, mutable_xzs, const_xzs, dim_x,
        order_right_left, phases, block_size, batch_size;
        workgroupsize = (block_SZ, 1), ndrange = (dim_x, dim_y)
        )
    if phase_B
        transform!(
            mod_4_identity!, mutable_phases, nothing; ndrange = dim_y
            )
    end
    return nothing

end

# CAUTION: Meta-programming is utilised to in order to avoid repetition.
for (safe_f_sym, unsafe_f_sym, right_left) in (
    (:mul_left!, :do_mul_left!, true), (:mul_right!, :do_mul_right!, false)
    )

#==============================================================================
RETURNS PAULI OPERATOR
==============================================================================#

# PauliOperator - PauliOperator
@eval @inline function $safe_f_sym(
    u::DevicePauliOperator, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    u.nqubits == v.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
    return $unsafe_f_sym(
        u, v;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
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

for (T_v_sym, v_tab_sym) in (
    (:DeviceTableau, :v), (:DeviceAbstractStabilizer, :(v.tab))
    )

# PauliOperator - Tableau/AbstractStabilizer[i]
@eval @inline function $safe_f_sym(
    u::DevicePauliOperator, v::$T_v_sym, i::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    1 <= i <= length($v_tab_sym.phases) || throw(BoundsError(THROW_BOUNDS))
    u.nqubits == $v_tab_sym.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
    return $unsafe_f_sym(
        u, v, i;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DevicePauliOperator, v::$T_v_sym, i::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    @inbounds device_mul!(
        u.phase, u.xz,
        (@view $v_tab_sym.phases[i]), (@view $v_tab_sym.xzs[:, i]);
        order_right_left = Val($right_left), phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Marks the end for (T_v_sym, v_tab_sym)
end

#==============================================================================
RETURNS TABLEAU / ABSTRACT STABILIZER
==============================================================================#

for (T_u_sym, u_tab_sym) in (
    (:DeviceTableau, :u), (:DeviceAbstractStabilizer, :(u.tab))
    )

# Tableau/AbstractStabilizer - PauliOperator
@eval @inline function $safe_f_sym(
    u::$T_u_sym, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    $u_tab_sym.nqubits == v.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
    return $unsafe_f_sym(
        u, v;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    device_mul!(
        $u_tab_sym.phases, $u_tab_sym.xzs, v.phase, v.xz;
        order_right_left = Val($right_left), phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau/AbstractStabilizer[i] - PauliOperator
@eval @inline function $safe_f_sym(
    u::$T_u_sym, i::Integer, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    1 <= i <= length($u_tab_sym.phases) || throw(BoundsError(THROW_BOUNDS))
    $u_tab_sym.nqubits == v.nqubits || throw(DimensionMismatch(THROW_NQUBITS))
    return $unsafe_f_sym(
        u, i, v;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, i::Integer, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    @inbounds device_mul!(
        (@view $u_tab_sym.phases[i]), (@view $u_tab_sym.xzs[:, i]),
        v.phase, v.xz;
        order_right_left = Val($right_left), phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# CAUTION: (Mixed)Destabilizer is handled separately.
# Tableau/AbstractStabilizer[i] - Self[j]
@eval @inline function $safe_f_sym(
    u::$T_u_sym, i::Integer, j::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    len = length($u_tab_sym.phases)
    1 <= i <= len && 1 <= j <= len || throw(BoundsError(THROW_BOUNDS))
    return $unsafe_f_sym(
        u, i, j;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, i::Integer, j::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    @inbounds device_mul!(
        (@view $u_tab_sym.phases[i]), (@view $u_tab_sym.xzs[:, i]),
        (@view $u_tab_sym.phases[j]), (@view $u_tab_sym.xzs[:, j]);
        order_right_left = Val($right_left), phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

for (T_v_sym, v_tab_sym) in (
    (:DeviceTableau, :v), (:DeviceAbstractStabilizer, :(v.tab))
    )

# Tableau/AbstractStabilizer - Tableau/AbstractStabilizer
@eval @inline function $safe_f_sym(
    u::$T_u_sym, v::$T_v_sym;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    length($u_tab_sym.phases) == length($v_tab_sym.phases) ||
        throw(DimensionMismatch(THROW_SIZE))
    $u_tab_sym.nqubits == $v_tab_sym.nqubits ||
        throw(DimensionMismatch(THROW_NQUBITS))
    return $unsafe_f_sym(
        u, v;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, v::$T_v_sym;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    device_mul!(
        $u_tab_sym.phases, $u_tab_sym.xzs, $v_tab_sym.phases, $v_tab_sym.xzs;
        order_right_left = Val($right_left), phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau/AbstractStabilizer - Tableau/AbstractStabilizer[i]
@eval @inline function $safe_f_sym(
    u::$T_u_sym, v::$T_v_sym, i::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    1 <= i <= length($v_tab_sym.phases) || throw(BoundsError(THROW_BOUNDS))
    $u_tab_sym.nqubits == $v_tab_sym.nqubits ||
        throw(DimensionMismatch(THROW_NQUBITS))
    return $unsafe_f_sym(
        u, v, i;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, v::$T_v_sym, i::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    @inbounds device_mul!(
        $u_tab_sym.phases, $u_tab_sym.xzs,
        (@view $v_tab_sym.phases[i]), (@view $v_tab_sym.xzs[:, i]);
        order_right_left = Val($right_left), phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau/AbstractStabilizer[i] - Tableau/AbstractStabilizer[j]
@eval @inline function $safe_f_sym(
    u::$T_u_sym, i::Integer, v::$T_v_sym, j::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    1 <= i <= length($u_tab_sym.phases) &&
        1 <= j <= length($v_tab_sym.phases) ||
            throw(BoundsError(THROW_BOUNDS))
    $u_tab_sym.nqubits == $v_tab_sym.nqubits ||
        throw(DimensionMismatch(THROW_NQUBITS))
    return $unsafe_f_sym(
        u, i, v, j;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, i::Integer, v::$T_v_sym, j::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    @inbounds device_mul!(
        (@view $u_tab_sym.phases[i]), (@view $u_tab_sym.xzs[:, i]),
        (@view $v_tab_sym.phases[j]), (@view $v_tab_sym.xzs[:, j]);
        order_right_left = Val($right_left), phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Marks the end for (T_v_sym, v_tab_sym)
end

# Marks the end for (T_u_sym, u_tab_sym)
end

#==============================================================================
RETURNS (MIXED) DESTABILIZER
==============================================================================#

# CAUTION: Requires special handling.
# (Mixed)Destabilizer[i] - Self[j]
@eval @inline function $safe_f_sym(
    u::DeviceUnionDestabilizer, i::Integer, j::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    len, n = length(u.tab.phases), u.tab.nqubits
    all(x -> 1 <= x <= len, (i, j, i + n, j + n)) ||
        throw(BoundsError(THROW_BOUNDS))
    return $unsafe_f_sym(
        u, i, j;
        phases = phases, block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceUnionDestabilizer, i::Integer, j::Integer;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    p, n, xzs = u.tab.phases, u.tab.nqubits, u.tab.xzs
    # Swapping the order of the indices is intentional.
    @inbounds device_mul!(
        (@view p[j]), (@view xzs[:, j]),
        (@view p[i]), (@view xzs[:, i]);
        order_right_left = Val($right_left), phases = Val(false),
        block_size = block_size, batch_size = batch_size
        )
    @inbounds device_mul!(
        (@view p[i + n]), (@view xzs[:, i + n]),
        (@view p[j + n]), (@view xzs[:, j + n]);
        order_right_left = Val($right_left), phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Marks the end for (safe_f_sym, unsafe_f_sym, right_left)
end
#=============================================================================#
