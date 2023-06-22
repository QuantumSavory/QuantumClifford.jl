module NoisyCircuits

#TODO permit the use of alternative RNGs

using QuantumClifford
using QuantumClifford: AbstractQCState, AbstractStabilizer, AbstractOperation, AbstractMeasurement, AbstractCliffordOperator, apply_single_x!, apply_single_y!, apply_single_z!, registered_statuses
import QuantumClifford: applywstatus!, affectedqubits

using Combinatorics: combinations
using Base.Cartesian

export AbstractOperation,
       UnbiasedUncorrelatedNoise,
       NoisyBellMeasurement,
       DecisionGate, ConditionalGate,
       affectedqubits, applynoise!,# TODO rename applyop to apply_mc
       applybranches, applynoise_branches,# TODO rename applybranches to apply_pe
       petrajectory, petrajectories,
       ConditionalGate, DecisionGate

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
struct DecisionGate <: AbstractOperation
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

"""Compute all possible new states after the application of the given operator. Reports the probability of each one of them. Deterministic, part of the Perturbative Expansion interface."""
function applybranches end

#TODO is the use of copy here necessary?
function applybranches(state, op; max_order=1)
    [(applywstatus!(copy(state),op)...,1,0)]
end

"""Compute all possible new states after the application of the given noise model. Reports the probability of each one of them. Deterministic, part of the Noise interface."""
function applynoise_branches end

function applynoise_branches(s::AbstractStabilizer,noise::UnbiasedUncorrelatedNoise,indices; max_order=1)
    n = nqubits(s)
    l = length(indices)
    infid = noise.errprobthird
    if l==0
        return [s,one(infid)]
    end
    error1 = 3*infid
    no_error1 = 1-error1
    no_error = no_error1^l
    results = [(copy(s),no_error,0)] # state, prob, order
    for order in 1:min(max_order,l)
        error_prob = no_error1^(l-order)*infid^order
        for error_indices in combinations(indices, order)
            _applynoise_branches_unbiased_uncorrelated(Val(order),s,error_indices,results,error_prob)
        end
    end
    results
end

@generated function _applynoise_branches_unbiased_uncorrelated(::Val{order},s,error_indices,results,error_prob) where {order}
    error_calls = Expr(:block)
    for i in 1:order
        call = quote (apply_single_x!,apply_single_y!,apply_single_z!)[$(Symbol(:i_,i))](new_state,error_indices[$i]) end
        push!(error_calls.args, call)
    end
    # n nested loops, one for each affected qubit, each loop dedicated to the 3 possible errors (X, Y, or Z)
    quote
        @nloops $order i d->1:3 begin
            new_state = copy(s)
            $error_calls
            push!(results,(new_state, error_prob, order))
        end
        results
    end
end

function applybranches(s::AbstractQCState, nop::NoiseOpAll; max_order=1)
    n = nqubits(s)
    return [(state, continue_stat, prob, order) for (state, prob, order) in applynoise_branches(s, nop.noise, 1:n, max_order=max_order)]
end

function applybranches(s::AbstractQCState, nop::NoiseOp; max_order=1)
    return [(state, continue_stat, prob, order) for (state, prob, order) in applynoise_branches(s, nop.noise, affectedqubits(nop), max_order=max_order)]
end

function applybranches(s::AbstractQCState, g::NoisyGate; max_order=1)
    news, _,_,_ = applybranches(s,g.gate,max_order=max_order)[1] # TODO this assumes only one always successful branch for the gate
    return [(state, continue_stat, prob, order) for (state, prob, order) in applynoise_branches(news, g.noise, affectedqubits(g), max_order=max_order)]
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

# TODO a lot of repetition with applywstatus!
function applybranches(s::AbstractQCState, m::BellMeasurement; max_order=1)
    n = nqubits(s)
    [(ns,iseven(r>>1 + m.parity) ? continue_stat : failure_stat, p,0)
     for (ns,r,p) in _applybranches_measurement([(s,0x0,1.0)],m.measurements,n)]
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

"""Run a perturbative expansion to a given order. Uses applybranches under the hood."""
function petrajectory(state, circuit; branch_weight=1.0, current_order=0, max_order=1)
    if size(circuit)[1] == 0
        return fill(zero(branch_weight), length(registered_statuses)-1)
    end
    next_op = circuit[1]
    rest_of_circuit = circuit[2:end]

    status_probs = fill(zero(branch_weight), length(registered_statuses)-1)

    for (i,(newstate, status, prob, order)) in enumerate(applybranches(state, next_op, max_order=max_order-current_order))
        if status==continue_stat # TODO is the copy below necessary?
            out_probs = petrajectory(copy(newstate), rest_of_circuit,
                branch_weight=branch_weight*prob, current_order=current_order+order, max_order=max_order)
            status_probs .+= out_probs
        else
            status_probs[status.status] += prob*branch_weight
        end
    end

    return status_probs
end

"""Run a perturbative expansion to a given order. This is the main public fuction for the perturbative expansion approach."""
function petrajectories(initialstate, circuit; branch_weight=1.0, max_order=1)
    status_probs = petrajectory(initialstate, circuit; branch_weight=branch_weight, current_order=0, max_order=max_order)
    Dict([CircuitStatus(i)=>status_probs[i] for i in eachindex(status_probs)])
end

function applynoise_branches(state::Register, noise, indices; max_order=1)
    [(Register(newstate,copy(state.bits)), prob, order)
     for (newstate, prob, order) in applynoise_branches(quantumstate(state), noise, indices; max_order=max_order)]
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

function QuantumClifford.apply!(state::Register, op::DecisionGate)
    decision = op.decisionfunction(state.bits)
    if !isnothing(decision)
        for i in 1:length(decision)
            apply!(state, op.gates[decision[i]])
        end
    end
    state
end


applybranches(s::Register, op::ConditionalGate; max_order=1) = [(applywstatus!(copy(s),op)...,1,0)]
applybranches(s::Register, op::DecisionGate; max_order=1) = [(applywstatus!(copy(s),op)...,1,0)]

function applybranches(s::Register, op::PauliMeasurement; max_order=1)
    stab = s.stab
    stab, anticom, r = project!(stab, op.pauli)
    new_branches = []
    if isnothing(r)
        s1 = s
        phases(stabilizerview(s1.stab))[anticom] = 0x00
        s1.bits[op.bit] = false
        s2 = copy(s)
        phases(stabilizerview(s2.stab))[anticom] = 0x02
        s2.bits[op.bit] = true
        push!(new_branches, (s1,continue_stat,1/2,0))
        push!(new_branches, (s2,continue_stat,1/2,0))
    else
        s.bits[op.bit] = r==0x02
        push!(new_branches, (s,continue_stat,1,0))
    end
    new_branches
end

#= TODO switch to sMX etc
function applybranches(state::Register, op::SparseMeasurement; max_order=1)
    n = nqubits(state.stab) # TODO implement actual sparse measurements
    p = zero(typeof(op.pauli), n)
    for (ii,i) in enumerate(op.indices)
        p[i] = op.pauli[ii]
    end
    dm = PauliMeasurement(p,op.bit)
    applybranches(state,dm, max_order=max_order)
end
=#

affectedqubits(m::NoisyBellMeasurement) = affectedqubits(m.meas)
affectedqubits(d::ConditionalGate) = union(affectedqubits(d.truegate), affectedqubits(d.falsegate))
affectedqubits(d::DecisionGate) = [(union(affectedqubits.(d.gates))...)...]
affectedqubits(v::VerifyOp) = v.indices

end
