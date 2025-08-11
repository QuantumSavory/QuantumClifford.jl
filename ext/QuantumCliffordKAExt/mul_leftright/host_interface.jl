
#=============================================================================#
# TODO: include the unsafe functions once the main package establishes them.
import QuantumClifford: mul_left!, mul_right!

# CAUTION: Meta-programming is utilised to in order to avoid repetition.
for (safe_f_sym, unsafe_f_sym, multiplication_order) in (
    (:mul_left!, :do_mul_left!, multiplication_order_left),
    (:mul_right!, :do_mul_right!, multiplication_order_right)
    )

#==============================================================================
RETURNS PAULI OPERATOR
==============================================================================#

# PauliOperator - PauliOperator
@eval @inline function $safe_f_sym(
    u::DevicePauliOperator, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        u.nqubits == v.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return $unsafe_f_sym(
        u, v;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DevicePauliOperator, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    device_mul!(
        u.phase, u.xz, v.phase, v.xz, $multiplication_order;
        phases = phases, primary_axis = primary_axis,
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
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        one(i) <= i <= length($v_tab_sym.phases) ||
            throw(BoundsError(THROW_BOUNDS))
        u.nqubits == $v_tab_sym.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return $unsafe_f_sym(
        u, v, i;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DevicePauliOperator, v::$T_v_sym, i::Integer;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @inbounds device_mul!(
        u.phase, u.xz,
        (@view $v_tab_sym.phases[i]), (@view $v_tab_sym.xzs[:, i]),
        $multiplication_order;
        phases = phases, primary_axis = primary_axis,
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
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        $u_tab_sym.nqubits == v.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return $unsafe_f_sym(
        u, v;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    device_mul!(
        $u_tab_sym.phases, $u_tab_sym.xzs, v.phase, v.xz,
        $multiplication_order;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau/AbstractStabilizer[i] - PauliOperator
@eval @inline function $safe_f_sym(
    u::$T_u_sym, i::Integer, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        one(i) <= i <= length($u_tab_sym.phases) ||
            throw(BoundsError(THROW_BOUNDS))
        $u_tab_sym.nqubits == v.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return $unsafe_f_sym(
        u, i, v;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, i::Integer, v::DevicePauliOperator;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @inbounds device_mul!(
        (@view $u_tab_sym.phases[i]), (@view $u_tab_sym.xzs[:, i]),
        v.phase, v.xz,
        $multiplication_order;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# CAUTION: (Mixed)Destabilizer is handled separately.
# Tableau/AbstractStabilizer[i] - Self[j]
@eval @inline function $safe_f_sym(
    u::$T_u_sym, i::Integer, j::Integer;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        len = length($u_tab_sym.phases)
        one(i) <= i <= len && one(j) <= j <= len ||
            throw(BoundsError(THROW_BOUNDS))
    end
    return $unsafe_f_sym(
        u, i, j;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, i::Integer, j::Integer;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @inbounds device_mul!(
        (@view $u_tab_sym.phases[i]), (@view $u_tab_sym.xzs[:, i]),
        (@view $u_tab_sym.phases[j]), (@view $u_tab_sym.xzs[:, j]),
        $multiplication_order;
        phases = phases, primary_axis = primary_axis,
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
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        length($u_tab_sym.phases) == length($v_tab_sym.phases) ||
            throw(DimensionMismatch(THROW_SIZE))
        $u_tab_sym.nqubits == $v_tab_sym.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return $unsafe_f_sym(
        u, v;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, v::$T_v_sym;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    device_mul!(
        $u_tab_sym.phases, $u_tab_sym.xzs, $v_tab_sym.phases, $v_tab_sym.xzs,
        $multiplication_order;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau/AbstractStabilizer - Tableau/AbstractStabilizer[i]
@eval @inline function $safe_f_sym(
    u::$T_u_sym, v::$T_v_sym, i::Integer;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        one(i) <= i <= length($v_tab_sym.phases) ||
            throw(BoundsError(THROW_BOUNDS))
        $u_tab_sym.nqubits == $v_tab_sym.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return $unsafe_f_sym(
        u, v, i;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, v::$T_v_sym, i::Integer;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @inbounds device_mul!(
        $u_tab_sym.phases, $u_tab_sym.xzs,
        (@view $v_tab_sym.phases[i]), (@view $v_tab_sym.xzs[:, i]),
        $multiplication_order;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau/AbstractStabilizer[i] - Tableau/AbstractStabilizer[j]
@eval @inline function $safe_f_sym(
    u::$T_u_sym, i::Integer, v::$T_v_sym, j::Integer;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        one(i) <= i <= length($u_tab_sym.phases) &&
            one(j) <= j <= length($v_tab_sym.phases) ||
                throw(BoundsError(THROW_BOUNDS))
        $u_tab_sym.nqubits == $v_tab_sym.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return $unsafe_f_sym(
        u, i, v, j;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::$T_u_sym, i::Integer, v::$T_v_sym, j::Integer;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @inbounds device_mul!(
        (@view $u_tab_sym.phases[i]), (@view $u_tab_sym.xzs[:, i]),
        (@view $v_tab_sym.phases[j]), (@view $v_tab_sym.xzs[:, j]),
        $multiplication_order;
        phases = phases, primary_axis = primary_axis,
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
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        len = length(u.tab.phases)
        n = len >> one(len)
        all(x -> one(x) <= x <= len, (i, j, i + n, j + n)) ||
            throw(BoundsError(THROW_BOUNDS))
    end
    return $unsafe_f_sym(
        u, i, j;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceUnionDestabilizer, i::Integer, j::Integer;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    p, xzs = u.tab.phases, u.tab.xzs
    n = length(p)
    n >>= one(n)
    # Swapping the order of the indices is intentional.
    @inbounds device_mul!(
        (@view p[j]), (@view xzs[:, j]),
        (@view p[i]), (@view xzs[:, i]),
        $multiplication_order;
        phases = Val(false), primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )
    @inbounds device_mul!(
        (@view p[i + n]), (@view xzs[:, i + n]),
        (@view p[j + n]), (@view xzs[:, j + n]),
        $multiplication_order;
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Marks the end for (safe_f_sym, unsafe_f_sym, multiplication_order)
end
#=============================================================================#
