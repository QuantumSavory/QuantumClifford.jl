
#=============================================================================#
# TODO: Include the unsafe functions once the main package establishes them.
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
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, v;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DevicePauliOperator, v::DevicePauliOperator;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, v;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

# PauliOperator - Tableau/AbstractStabilizer[i]
@eval @inline function $safe_f_sym(
    u::DevicePauliOperator, v::DeviceUnionTableau, i::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, v, i;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DevicePauliOperator, v::DeviceUnionTableau, i::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, v, i;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

#==============================================================================
RETURNS TABLEAU / ABSTRACT STABILIZER
==============================================================================#

# Tableau/AbstractStabilizer - PauliOperator
@eval @inline function $safe_f_sym(
    u::DeviceUnionTableau, v::DevicePauliOperator;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, v;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceUnionTableau, v::DevicePauliOperator;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, v;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

# Tableau/AbstractStabilizer - Tableau/AbstractStabilizer[i]
@eval @inline function $safe_f_sym(
    u::DeviceUnionTableau, v::DeviceUnionTableau, i::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, v, i;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceUnionTableau, v::DeviceUnionTableau, i::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, v, i;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

# Tableau - Tableau/AbstractStabilizer
@eval @inline function $safe_f_sym(
    u::DeviceTableau, v::DeviceUnionTableau;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, v;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceTableau, v::DeviceUnionTableau;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, v;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

# Tableau[i] - PauliOperator
@eval @inline function $safe_f_sym(
    u::DeviceTableau, i::Integer, v::DevicePauliOperator;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, i, v;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceTableau, i::Integer, v::DevicePauliOperator;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, i, v;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

# Tableau[i] - Tableau/AbstractStabilizer[j]
@eval @inline function $safe_f_sym(
    u::DeviceTableau, i::Integer, v::DeviceUnionTableau, j::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, i, v, j;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceTableau, i::Integer, v::DeviceUnionTableau, j::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, i, v, j;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

# CAUTION: (Mixed)Destabilizer is handled separately.
# Tableau/AbstractStabilizer[i] - Self[j]
@eval @inline function $safe_f_sym(
    u::DeviceUnionTableau, i::Integer, j::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, i, j;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceUnionTableau, i::Integer, j::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, i, j;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

#==============================================================================
RETURNS (MIXED) DESTABILIZER
==============================================================================#

# CAUTION: Requires special handling.
# (Mixed)Destabilizer[i] - Self[j]
@eval @inline function $safe_f_sym(
    u::DeviceUnionDestabilizer, i::Integer, j::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return mul!(
        u, i, j;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

@eval @inline function $unsafe_f_sym(
    u::DeviceUnionDestabilizer, i::Integer, j::Integer;
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Val{phases_B} = Val(default_phases),
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    ) where {phases_B}

    return do_mul!(
        u, i, j;
        multiplication_order = $multiplication_order,
        primary_axis = primary_axis, phases = phases_B,
        block_size = block_size, batch_size = batch_size
        )

end

# Marks the end for (safe_f_sym, unsafe_f_sym, multiplication_order)
end
#=============================================================================#
