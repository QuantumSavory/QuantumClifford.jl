
#=============================================================================#
# PauliOperator
function adapt_structure(
    AT::Type{T}, pauli::PauliOperator
    ) where {T <: AbstractArray}

    return PauliOperator(
        adapt(AT, change_type(UInt8, pauli.phase)),
        pauli.nqubits,
        adapt(AT, pauli.xz)
        )

end

function adapt_structure(
    AT::Type{T}, pauli::PauliOperator
    ) where {T <: AbstractGPUArray}

    return PauliOperator(
        adapt(AT, change_type(DeviceUnsigned, pauli.phase)),
        pauli.nqubits,
        adapt(AT, pauli.xz)
        )

end

# Tableau
function adapt_structure(
    AT::Type{T}, tab::Tableau
    ) where {T <: AbstractArray}

    return Tableau(
        adapt(AT, change_type(UInt8, tab.phases)),
        tab.nqubits,
        adapt(AT, tab.xzs)
        )

end

function adapt_structure(
    AT::Type{T}, tab::Tableau
    ) where {T <: AbstractGPUArray}

    return Tableau(
        adapt(AT, change_type(DeviceUnsigned, tab.phases)),
        tab.nqubits,
        adapt(AT, tab.xzs)
        )

end

# Stabilizer
function adapt_structure(
    AT::Type{T}, state::Stabilizer
    ) where {T <: AbstractArray}

    return Stabilizer(adapt(AT, state.tab))

end

# Destabilizer
function adapt_structure(
    AT::Type{T}, state::Destabilizer
    ) where {T <: AbstractArray}

    return Destabilizer(adapt(AT, state.tab))

end

# MixedStabilizer
function adapt_structure(
    AT::Type{T}, state::MixedStabilizer
    ) where {T <: AbstractArray}

    return MixedStabilizer(adapt(AT, state.tab), state.rank)

end

# MixedDestabilizer
function adapt_structure(
    AT::Type{T}, state::MixedDestabilizer
    ) where {T <: AbstractArray}

    return MixedDestabilizer(adapt(AT, state.tab), state.rank)

end
#=============================================================================#
