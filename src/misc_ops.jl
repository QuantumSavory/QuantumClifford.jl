import QuantumInterface: nsubsystems

"""A Stabilizer measurement on the entirety of the quantum register.

`projectrand!(state, pauli)` and `apply!(state, PauliMeasurement(pauli))` give the same (possibly non-deterministic) result.
Particularly useful when acting on [`Register`](@ref).

See also: [`apply!`](@ref), [`projectrand!`](@ref)."""
struct PauliMeasurement{
    P <: AbstractArray{<: Unsigned, 0}, XZ <: AbstractVector{<: Unsigned}
} <: AbstractMeasurement
    pauli::PauliOperator{P,XZ}
    bit::Int
end

PauliMeasurement(pauli) = PauliMeasurement(pauli,0)
PauliMeasurement(pauli,::Nothing) = PauliMeasurement(pauli,0)

function apply!(state::AbstractStabilizer, m::PauliMeasurement)
    projectrand!(state,m.pauli)
    state
end

function apply!(state::MixedDestabilizer, indices::Base.AbstractVecOrTuple, operation::Type{<:AbstractSymbolicOperator})
    apply!(state, operation(indices...))
end

"""A Clifford gate, applying the given `cliff` operator to the qubits at the selected `indices`.

`apply!(state, cliff, indices)` and `apply!(state, SparseGate(cliff, indices))` give the same result."""
struct SparseGate{T<:Tableau} <: AbstractCliffordOperator # TODO simplify type parameters (remove nesting)
    cliff::CliffordOperator{T}
    indices::Vector{Int}
    function SparseGate(cliff::CliffordOperator{T}, indices::Vector{Int}) where T<:Tableau
        if length(indices) != nqubits(cliff)
            throw(ArgumentError("The number of target qubits (indices) must match the qubit count in the CliffordOperator."))
        end
        new{T}(cliff, indices)
    end
end

SparseGate(c,t::Tuple) = SparseGate(c,collect(t))

function apply!(state::AbstractStabilizer, g::SparseGate; kwargs...)
    m = maximum(g.indices)
    if m > nqubits(state)
        throw(ArgumentError(lazy"SparseGate was attempted on invalid qubit index $(m) when the state contains only $(nqubits(state)) qubits."))
    end
    apply!(state, g.cliff, g.indices; kwargs...)
end

function LinearAlgebra.inv(g::SparseGate; phases=true)
  return SparseGate(inv(g.cliff;phases=phases), g.indices)
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

"""A Bell measurement in which each of the measured qubits has a chance to have flipped."""
struct NoisyBellMeasurement{T} <: AbstractOperation
    meas::AbstractOperation
    flipprob::T
end
NoisyBellMeasurement(p,i,fp) = NoisyBellMeasurement(BellMeasurement(p,i),fp)

function applywstatus!(s::AbstractQCState, m::NoisyBellMeasurement)
    state, status = applywstatus!(s,m.meas)
    nqubits = length(affectedqubits(m))
    errprob = (1-(1-2*m.flipprob)^nqubits)/2 # probability of odd number of flips
    if rand()<errprob
        return state, status==continue_stat ? failure_stat : continue_stat
    else
        return state, status
    end
end

function applybranches(s::AbstractQCState, m::BellMeasurement; max_order=1)
    n = nqubits(s)
    [(ns,iseven(r>>1 + m.parity) ? continue_stat : failure_stat, p,0)
     for (ns,r,p) in _applybranches_measurement([(s,0x0,1.0)],m.measurements,n)]
end

function applybranches(s::AbstractQCState, m::NoisyBellMeasurement; max_order=1)
    measurement_branches = applybranches(s, m.meas, max_order=max_order)
    if max_order==0
        return measurement_branches
    else
        new_branches = []
        nqubits = length(affectedqubits(m))
        p = (1-2m.flipprob)^nqubits
        errprob = 1//2*(1-p)
        sucprob = 1//2*(1+p)
        for (mstate, success, mprob, morder) in measurement_branches
            push!(new_branches, (mstate, success, mprob*sucprob, morder))
            push!(new_branches, (mstate, success==continue_stat ? failure_stat : continue_stat, mprob*errprob, morder+1))
        end
        return new_branches
    end
end

# TODO XXX THIS IS PARTICULARLY INEFFICIENT recurrent implementation
function _applybranches_measurement(branches, measurements, n)
    if length(measurements) == 0
        return branches
    end

    new_branches = []
    pauli = measurements[1]
    otherpaulis = measurements[2:end]

    for (s,r0,p) in branches
        _,anticom,r = project!(quantumstate(s),pauli)
        if isnothing(r) # TODO anticom could be zero if there was a rank change
            s1 = s
            s2 = copy(s)
            r1 = phases(stabilizerview(s1))[anticom] = 0x00
            r2 = phases(stabilizerview(s2))[anticom] = 0x02
            push!(new_branches, (s1,r0+r1,p/2))
            push!(new_branches, (s2,r0+r2,p/2))
        else
            push!(new_branches, (s,r0+r,p))
        end
    end

    return _applybranches_measurement(new_branches, otherpaulis, n)
end


"""A "probe" to verify that the state of the qubits corresponds to a desired `good_state`, e.g. at the end of the execution of a circuit."""
struct VerifyOp <: AbstractOperation
    good_state::Stabilizer
    indices::AbstractVector{Int}
    function VerifyOp(s, indices)
        gs = canonicalize_rref!(copy(stabilizerview(s)))[1]
        r,n = size(gs)
        if(r!=n)
            throw(ArgumentError(lazy"""The argument you have provided for good_state is not a logical state within the codespace. Expected a pure $n - qubit stabilizer state (i.e. $n independent stabilizer generators on $n qubits), but good_state has only $r independent stabilizer generators."""))
        end
        return new(gs, indices)
    end
end

# TODO this one needs more testing
function applywstatus!(s::AbstractQCState, v::VerifyOp) # XXX It assumes the other qubits are measured or traced out
    # TODO QuantumClifford should implement some submatrix comparison
    sv = traceout!(copy(quantumstate(s)), setdiff(1:nqubits(s),v.indices))
    sv = stabilizerview(sv)
    canonicalize_rref!(sv, v.indices)
    sv = tab(sv)
    good_state = tab(v.good_state)
    
    for i in eachindex(good_state)
        (sv.phases[end-i+1]==good_state.phases[end-i+1]) || return s, false_success_stat
        for (j,q) in zip(eachindex(good_state),v.indices)
            (sv[end-i+1,q]==good_state[end-i+1,j]) || return s, false_success_stat
        end
    end
    return s, true_success_stat
end

operatordeterminism(::Type{VerifyOp}) = DeterministicOperatorTrait()

"""Applies an XOR gate to classical bits. Currently only implemented for functionality with pauli frames."""
struct ClassicalXOR{N} <: AbstractOperation
    "The indices of the classical bits to be xor-ed"
    bits::NTuple{N,Int}
    "The index of the classical bit that will store the results"
    store::Int
end

ClassicalXOR(bits,store) = ClassicalXOR{length(bits)}(tuple(bits...),store)

nsubsystems(state::MixedDestabilizer) = nqubits(state)
