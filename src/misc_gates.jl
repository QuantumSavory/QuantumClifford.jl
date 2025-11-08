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

affectedqubits(d::ConditionalGate) = union(affectedqubits(d.truegate), affectedqubits(d.falsegate))
affectedqubits(d::IndexedDecisionGate) = [(union(affectedqubits.(d.gates))...)...]