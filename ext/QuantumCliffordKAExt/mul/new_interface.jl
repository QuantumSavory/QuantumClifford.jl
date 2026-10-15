
#=============================================================================#
# TODO: Import the functions once the main package establishes them.

#==============================================================================
RETURNS PAULI OPERATOR
==============================================================================#

# PauliOperator - PauliOperator
@inline function mul!(
    u::DevicePauliOperator, v::DevicePauliOperator;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        u.nqubits == v.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return do_mul!(
        u, v;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DevicePauliOperator, v::DevicePauliOperator;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    device_mul!(
        u.phase, u.xz,
        v.phase, v.xz;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# PauliOperator - Tableau/AbstractStabilizer[i]
@inline function mul!(
    u::DevicePauliOperator, v::DeviceUnionTableau, i::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        v_tab = tab(v)
        one(i) <= i <= length(v_tab.phases) ||
            throw(BoundsError(THROW_BOUNDS))
        u.nqubits == v_tab.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return do_mul!(
        u, v, i;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DevicePauliOperator, v::DeviceUnionTableau, i::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    v_tab = tab(v)
    @inbounds device_mul!(
        u.phase, u.xz,
        (@view v_tab.phases[i]), (@view v_tab.xzs[:, i]);
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

#==============================================================================
RETURNS TABLEAU / ABSTRACT STABILIZER
==============================================================================#

# Tableau/AbstractStabilizer - PauliOperator
@inline function mul!(
    u::DeviceUnionTableau, v::DevicePauliOperator;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        u_tab = tab(u)
        u_tab.nqubits == v.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return do_mul!(
        u, v;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DeviceUnionTableau, v::DevicePauliOperator;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    u_tab = tab(u)
    device_mul!(
        u_tab.phases, u_tab.xzs,
        v.phase, v.xz;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau/AbstractStabilizer - Tableau/AbstractStabilizer[i]
@inline function mul!(
    u::DeviceUnionTableau, v::DeviceUnionTableau, i::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        u_tab, v_tab = tab(u), tab(v)
        one(i) <= i <= length(v_tab.phases) ||
            throw(BoundsError(THROW_BOUNDS))
        u_tab.nqubits == v_tab.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return do_mul!(
        u, v, i;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DeviceUnionTableau, v::DeviceUnionTableau, i::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    u_tab, v_tab = tab(u), tab(v)
    @inbounds device_mul!(
        u_tab.phases, u_tab.xzs,
        (@view v_tab.phases[i]), (@view v_tab.xzs[:, i]);
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau - Tableau/AbstractStabilizer
@inline function mul!(
    u::DeviceTableau, v::DeviceUnionTableau;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        v_tab = tab(v)
        length(u.phases) == length(v_tab.phases) ||
            throw(DimensionMismatch(THROW_SIZE))
        u.nqubits == v_tab.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return do_mul!(
        u, v;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DeviceTableau, v::DeviceUnionTableau;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    v_tab = tab(v)
    device_mul!(
        u.phases, u.xzs,
        v_tab.phases, v_tab.xzs;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau[i] - PauliOperator
@inline function mul!(
    u::DeviceTableau, i::Integer, v::DevicePauliOperator;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        one(i) <= i <= length(u.phases) ||
            throw(BoundsError(THROW_BOUNDS))
        u.nqubits == v.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return do_mul!(
        u, i, v;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DeviceTableau, i::Integer, v::DevicePauliOperator;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @inbounds device_mul!(
        (@view u.phases[i]), (@view u.xzs[:, i]),
        v.phase, v.xz;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# Tableau[i] - Tableau/AbstractStabilizer[j]
@inline function mul!(
    u::DeviceTableau, i::Integer, v::DeviceUnionTableau, j::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        v_tab = tab(v)
        one(i) <= i <= length(u.phases) &&
            one(j) <= j <= length(v_tab.phases) ||
                throw(BoundsError(THROW_BOUNDS))
        u.nqubits == v_tab.nqubits ||
            throw(DimensionMismatch(THROW_NQUBITS))
    end
    return do_mul!(
        u, i, v, j;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DeviceTableau, i::Integer, v::DeviceUnionTableau, j::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    v_tab = tab(v)
    @inbounds device_mul!(
        (@view u.phases[i]), (@view u.xzs[:, i]),
        (@view v_tab.phases[j]), (@view v_tab.xzs[:, j]);
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

# CAUTION: (Mixed)Destabilizer is handled separately.
# Tableau/AbstractStabilizer[i] - Self[j]
@inline function mul!(
    u::DeviceUnionTableau, i::Integer, j::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        u_tab = tab(u)
        len = length(u_tab.phases)
        one(i) <= i <= len && one(j) <= j <= len ||
            throw(BoundsError(THROW_BOUNDS))
    end
    return do_mul!(
        u, i, j;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DeviceUnionTableau, i::Integer, j::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    u_tab = tab(u)
    @inbounds device_mul!(
        (@view u_tab.phases[i]), (@view u_tab.xzs[:, i]),
        (@view u_tab.phases[j]), (@view u_tab.xzs[:, j]);
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end

#==============================================================================
RETURNS (MIXED) DESTABILIZER
==============================================================================#

# CAUTION: Requires special handling.
# (Mixed)Destabilizer[i] - Self[j]
@inline function mul!(
    u::DeviceUnionDestabilizer, i::Integer, j::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        len = length(u.tab.phases) >> 1
        one(i) <= i <= len && one(j) <= j <= len ||
            throw(BoundsError(THROW_BOUNDS))
    end
    return do_mul!(
        u, i, j;
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_mul!(
    u::DeviceUnionDestabilizer, i::Integer, j::Integer;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    len = length(u.tab.phases) >> 1
    # Swapping the order of the indices is intentional.
    @inbounds device_mul!(
        (@view u.tab.phases[j]), (@view u.tab.xzs[:, j]),
        (@view u.tab.phases[i]), (@view u.tab.xzs[:, i]);
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = false,
        block_size = block_size, batch_size = batch_size
        )
    @inbounds device_mul!(
        (@view u.tab.phases[i + len]), (@view u.tab.xzs[:, i + len]),
        (@view u.tab.phases[j + len]), (@view u.tab.xzs[:, j + len]);
        multiplication_order = multiplication_order,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )
    return u

end
#=============================================================================#
