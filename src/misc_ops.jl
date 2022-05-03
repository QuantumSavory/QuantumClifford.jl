"""A Stabilizer measurement on the entirety of the quantum register.

`projectrand!(state, pauli)` and `apply!(state, PauliMeasurement(pauli))` give the same (possibly non-deterministic) result.
Particularly useful when acting on [`Register`](@ref).

See also: [`apply!`](@ref), [`projectrand!`](@ref)."""
struct PauliMeasurement{Tz<:AbstractArray{UInt8,0}, Tv<:AbstractVector{<:Unsigned}, T<:Union{Int,Nothing}} <: AbstractMeasurement
    pauli::PauliOperator{Tz,Tv}
    storagebit::T
end

PauliMeasurement(pauli) = PauliMeasurement(pauli,nothing)

function apply!(state::AbstractStabilizer, m::PauliMeasurement)
    projectrand!(state,m.pauli)
    state
end

"""A Clifford gate, applying the given `cliff` operator to the qubits at the selected `indices`.

`apply!(state, cliff, indices)` and `apply!(state, SparseGate(cliff, indices))` give the same result."""
struct SparseGate{Tzv<:AbstractVector{UInt8},Tm<:AbstractMatrix{<:Unsigned}} <: AbstractCliffordOperator
    cliff::CliffordOperator{Tzv,Tm}
    indices::Vector{Int}
end

function apply!(state::AbstractStabilizer, g::SparseGate; kwargs...)
    apply!(state, g.cliff, g.indices; kwargs...)
end

"""Reset the specified qubits to the given state."""
struct Reset{Tzv<:AbstractVector{UInt8},Tm<:AbstractMatrix{<:Unsigned}} <: AbstractOperation 
    resetto::Stabilizer{Tzv, Tm}
    indices::Vector{Int}
end

function apply!(state::AbstractStabilizer, reset::Reset)
    reset_qubits!(state, reset.resetto, reset.indices)
    return s
end