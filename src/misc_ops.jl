"""A Stabilizer measurement on the entirety of the quantum register.

`projectrand!(state, pauli)` and `apply!(state, PauliMeasurement(pauli))` give the same (possibly non-deterministic) result.
Particularly useful when acting on [`Register`](@ref).

See also: [`apply!`](@ref), [`projectrand!`](@ref)."""
struct PauliMeasurement{Tz<:AbstractArray{UInt8,0}, Tv<:AbstractVector{<:Unsigned}} <: AbstractMeasurement
    pauli::PauliOperator{Tz,Tv}
    bit::Int
end

PauliMeasurement(pauli) = PauliMeasurement(pauli,0)
PauliMeasurement(pauli,::Nothing) = PauliMeasurement(pauli,0)

function apply!(state::AbstractStabilizer, m::PauliMeasurement)
    projectrand!(state,m.pauli)
    state
end

"""A Clifford gate, applying the given `cliff` operator to the qubits at the selected `indices`.

`apply!(state, cliff, indices)` and `apply!(state, SparseGate(cliff, indices))` give the same result."""
struct SparseGate{T<:Tableau} <: AbstractCliffordOperator # TODO simplify type parameters (remove nesting)
    cliff::CliffordOperator{T}
    indices::Vector{Int}
end

SparseGate(c,t::Tuple) = SparseGate(c,collect(t))

function apply!(state::AbstractStabilizer, g::SparseGate; kwargs...)
    apply!(state, g.cliff, g.indices; kwargs...)
end

"""Reset the specified qubits to the given state.

Be careful, this operation implies first tracing out the qubits, which can lead to mixed states
if these qubits were entangled with the rest of the system.

See also: [`sMRZ`](@ref)"""
struct Reset{T<:Tableau} <: AbstractOperation # TODO simplify type parameters (remove nesting)
    resetto::Stabilizer{T}
    indices::Vector{Int}
end

function apply!(state::AbstractStabilizer, reset::Reset)
    reset_qubits!(state, reset.resetto, reset.indices)
    return state
end

"""A Bell measurement performing the correlation measurement corresponding to the given `pauli` projections on the qubits at the selected indices."""
struct BellMeasurement <: AbstractOperation
    measurements::Vector{Union{sMX,sMY,sMZ}}
    parity::Bool
end
BellMeasurement(meas) = BellMeasurement(meas,false)

# TODO this seems unnecessarily complicated
function applywstatus!(s::AbstractQCState, m::BellMeasurement)
    res = 0x00
    for meas in m.measurements
        s,r = projectrand!(s,meas)
        res ⊻= r
    end
    if iseven(res>>1+m.parity)
        return s, continue_stat
    else
        return s, failure_stat
    end
end

"""A "probe" to verify that the state of the qubits corresponds to a desired `good_state`, e.g. at the end of the execution of a circuit."""
struct VerifyOp <: AbstractOperation
    good_state::Stabilizer
    indices::AbstractVector{Int}
    VerifyOp(s,indices) = new(canonicalize_rref!(copy(stabilizerview(s)))[1],indices)
end

# TODO this one needs more testing
function applywstatus!(s::AbstractQCState, v::VerifyOp) # XXX It assumes the other qubits are measured or traced out
    # TODO QuantumClifford should implement some submatrix comparison
    canonicalize_rref!(quantumstate(s),v.indices) # Document why rref is used
    sv = tab(s)
    good_state = tab(v.good_state)
    for i in eachindex(good_state)
        (sv.phases[end-i+1]==good_state.phases[end-i+1]) || return s, false_success_stat
        for (j,q) in zip(eachindex(good_state),v.indices)
            (sv[end-i+1,q]==good_state[end-i+1,j]) || return s, false_success_stat
        end
    end
    return s, true_success_stat
end
