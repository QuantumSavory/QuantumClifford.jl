module NoisyCircuits

#TODO permit the use of alternative RNGs

using QuantumClifford
using QuantumClifford: AbstractQCState, AbstractStabilizer, AbstractOperation, AbstractMeasurement, AbstractCliffordOperator, apply_single_x!, apply_single_y!, apply_single_z!, registered_statuses
import QuantumClifford: applywstatus!, affectedqubits, applybranches, applynoise_branches

using Combinatorics: combinations
using Base.Cartesian

export NoisyBellMeasurement,
       IndexedDecisionGate, ConditionalGate

#TODO all these structs should use specified types
#TODO all of the methods should better specified type signatures
#TODO measure allocation in various apply* methods and verify it is not superfluous

"""A Bell measurement in which each of the measured qubits has a chance to have flipped."""
struct NoisyBellMeasurement{T} <: AbstractOperation
    meas::AbstractOperation
    flipprob::T
end
NoisyBellMeasurement(p,i,fp) = NoisyBellMeasurement(BellMeasurement(p,i),fp)

"""A conditional gate that either performs `truegate` or `falsegate`, depending on the value of `controlbit`."""
struct ConditionalGate <: AbstractOperation
    truegate::AbstractOperation
    falsegate::AbstractOperation
    controlbit::Int
end

"""A conditional gate that performs one of the `gates`, depending on the output of `decisionfunction` applied to the entire classical bit register."""
struct IndexedDecisionGate <: AbstractOperation
    gates::AbstractVector{AbstractOperation}
    decisionfunction
end

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

# TODO tests for this
function QuantumClifford.apply!(state::Register, op::ConditionalGate)
    if state.bits[op.controlbit]
        apply!(state, op.truegate)
    else
        apply!(state, op.falsegate)
    end
    return state
end

function QuantumClifford.apply!(state::Register, op::IndexedDecisionGate)
    decision = op.decisionfunction(state.bits)
    if !isnothing(decision)
        for i in 1:length(decision)
            apply!(state, op.gates[decision[i]])
        end
    end
    state
end

applybranches(s::Register, op::ConditionalGate; max_order=1) = [(applywstatus!(copy(s),op)...,1,0)]
applybranches(s::Register, op::IndexedDecisionGate; max_order=1) = [(applywstatus!(copy(s),op)...,1,0)]

affectedqubits(m::NoisyBellMeasurement) = affectedqubits(m.meas)
affectedqubits(d::ConditionalGate) = union(affectedqubits(d.truegate), affectedqubits(d.falsegate))
affectedqubits(d::IndexedDecisionGate) = [(union(affectedqubits.(d.gates))...)...]

end
