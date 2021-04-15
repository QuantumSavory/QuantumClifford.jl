"""
A module for simulating noisy Clifford circuits.
"""
module NoisyCircuits

#TODO permit the use of alternative RNGs

using QuantumClifford
using QuantumClifford: AbstractStabilizer, AbstractCliffordOperator

using StatsBase: countmap
using Combinatorics: combinations

export AbstractOperation,
       UnbiasedUncorrelatedNoise, NoiseOp, NoiseOpAll, VerifyOp,
       SparseGate, DenseGate, NoisyGate,
       BellMeasurement, NoisyBellMeasurement,
       DenseMeasurement,
       DecisionGate, ConditionalGate,
       affectedqubits, applyop!, applynoise!,
       applyop_branches, applynoise_branches,
       mctrajectory!, mctrajectories,
       petrajectory, petrajectories,
       Register, Measurement, ConditionalGate, DecisionGate

abstract type AbstractOperation end

abstract type AbstractNoise end

#TODO all these structs should use specified types
#TODO all of the methods should better specified type signatures

"""Depolarization noise model with total probability of error `3*errprobthird`."""
struct UnbiasedUncorrelatedNoise{T} <: AbstractNoise
    errprobthird::T
end

"""An operator that applies the given `noise` model to the qubits at the selected `indices`."""
struct NoiseOp <: AbstractOperation
    noise::AbstractNoise
    indices::AbstractVector{Int}
end

"""An operator that applies the given `noise` model to all qubits."""
struct NoiseOpAll <: AbstractOperation
    noise::AbstractNoise
end

"""A Clifford gate, applying the given `cliff` operator to the qubits at the selected `indices`."""
struct SparseGate <: AbstractOperation
    cliff::AbstractCliffordOperator
    indices::AbstractVector{Int}
end

"""A Clifford gate, applying the given `cliff` operator."""
struct DenseGate <: AbstractOperation
    cliff::AbstractCliffordOperator
end

"""A gate consisting of the given `noise` applied after the given perfect Clifford `gate`."""
struct NoisyGate <: AbstractOperation
    gate::AbstractOperation
    noise::AbstractNoise
end

"""A Bell measurement performing the correlation measurement corresponding to the given `pauli` projections on the qubits at the selected indices."""
struct BellMeasurement <: AbstractOperation
    pauli::AbstractVector{PauliOperator}
    indices::AbstractVector{Int}
end

"""A Bell measurement in which each of the measured qubits has a chance to have flipped."""
struct NoisyBellMeasurement{T} <: AbstractOperation
    meas::AbstractOperation
    flipprob::T
end

"""Performing a Bell measurement followed by resetting the measured qubits to the given state `resetto`."""
struct BellMeasurementAndReset <: AbstractOperation # TODO delete this and just make a standalone reset 
    meas::AbstractOperation
    resetto::Stabilizer
end

"""A "probe" to verify that the state of the qubits corresponds to a desired `good_state`, e.g. at the end of the execution of a circuit."""
struct VerifyOp <: AbstractOperation
    good_state::Stabilizer
    indices::AbstractVector{Int}
    VerifyOp(s,indices) = new(canonicalize_rref!(copy(s))[1],indices)
end

"""A Stabilizer measurement on the """
struct DenseMeasurement <: AbstractOperation
    pauli::PauliOperator
    storagebit::Int
end

struct ConditionalGate <: AbstractOperation
    truegate::AbstractOperation
    falsegate::AbstractOperation
    controlbit::Int
end

struct DecisionGate <: AbstractOperation
    gates::AbstractVector{AbstractOperation}
    decisionfunction
end

"""A list of default statuses returned by `applyop!`."""
const statuses = [:continue, :detected_failure, :undetected_failure, :true_success]

"""A method giving the qubits acted upon by a given operation. Part of the Noise interface."""
function affectedqubits end
affectedqubits(g::NoisyGate) = affectedqubits(g.gate)
affectedqubits(g::SparseGate) = g.indices
affectedqubits(g::DenseGate) = 1:nqubits(g.cliff)
affectedqubits(m::BellMeasurement) = m.indices
affectedqubits(m::BellMeasurementAndReset) = affectedqubits(m.meas)
affectedqubits(m::NoisyBellMeasurement) = affectedqubits(m.meas)
affectedqubits(n::NoiseOp) = n.indices
affectedqubits(v::VerifyOp) = v.indices
affectedqubits(g::DenseMeasurement) = 1:length(g.pauli)
affectedqubits(d::ConditionalGate) = union(affectedqubits(d.truegate), affectedqubits(d.falsegate))
affectedqubits(d::DecisionGate) = [(union(affectedqubits.(d.gates))...)...]

"""A method modifying a given state by applying the given operation. Non-deterministic, part of the Monte Carlo interface."""
function applyop! end

function applyop!(s::Stabilizer, g::NoisyGate)
    s = applynoise!(
            applyop!(s,g.gate)[1],
            g.noise,
            affectedqubits(g.gate)),
    return s, :continue
end

applyop!(s::Stabilizer, g::SparseGate) = (apply!(s,g.cliff,affectedqubits(g)), :continue)

applyop!(s::Stabilizer, g::DenseGate) = (apply!(s,g.cliff), :continue)

function applyop!(s::Stabilizer, m::NoisyBellMeasurement)
    state, status = applyop!(s,m.meas)
    nqubits = length(affectedqubits(m))
    errprob = (1-(1-2*m.flipprob)^nqubits)/2 # probability of odd number of flips
    if rand()<errprob
        return state, status==:continue ? :detected_failure : :continue
    else
        return state, status
    end
end

# TODO this seems unnecessarily complicated
function applyop!(s::Stabilizer, m::BellMeasurement)
    n = nqubits(s)
    indices = affectedqubits(m)
    res = 0x00
    for (pauli, index) in zip(m.pauli,affectedqubits(m))
        if pauli == X # TODO this is not an elegant way to choose between X and Z coincidence measurements
            op = single_x(n,index) # TODO this is pretty terribly inefficient... use some sparse check
        elseif pauli == Z
            op = single_z(n,index)
        elseif pauli == Y
            op = single_y(n,index)
        else
            op = -single_y(n, index)
        end
         # TODO permit Y operators and permit negative operators
        s,anticom,r = project!(s,op)
        if isnothing(r)
            if rand()>0.5 # TODO this seems stupid, float not necessary
                r = s.phases[anticom] = 0x00
            else
                r = s.phases[anticom] = 0x02
            end
        end
        res ⊻= r
    end
    if res==0x0
        return s, :continue
    else
        return s, :detected_failure
    end
end

function applyop!(s::Stabilizer, mr::BellMeasurementAndReset)
    s,res = applyop!(s,mr.meas)
    if !res
        return s,res
    else
        # TODO is the traceout necessary given that we just performed measurements?
        traceout!(s,mr.meas.indices)# TODO it seems like a bad idea not to keep track of the rank here
        for (ii,i) in enumerate(affectedqubits(mr))
            for j in [1,2]
                s[end-j+1,i] = mr.resetto[j,ii]
            end
        end
        return s,:continue
    end
end

"""A method modifying a given state by applying the corresponding noise model. Non-deterministic, part of the Noise interface."""
function applynoise! end

function applynoise!(s::Stabilizer,noise::UnbiasedUncorrelatedNoise,indices::AbstractVector{Int})
    n = nqubits(s)
    infid = noise.errprobthird
    for i in indices
        r = rand()
        if r<infid
            apply_single_x!(s,i)
        end
        if infid<=r<2infid
            apply_single_z!(s,i)
        end
        if 2infid<=r<3infid
            apply_single_y!(s,i)
        end
    end
    s
end

function applyop!(s::Stabilizer, mr::NoiseOpAll)
    n = nqubits(s)
    return applynoise!(s, mr.noise, 1:n), :continue
end

function applyop!(s::Stabilizer, mr::NoiseOp)
    return applynoise!(s, mr.noise, affectedqubits(mr)), :continue
end

# TODO this one needs more testing
function applyop!(s::Stabilizer, v::VerifyOp) # XXX It assumes the other qubits are measured or traced out
    # TODO QuantumClifford should implement some submatrix comparison
    s, _ = canonicalize_rref!(s,v.indices) # Document why rref is used
    for i in eachindex(v.good_state)
        (s.phases[end-i+1]==v.good_state.phases[end-i+1]) || return s, :undetected_failure
        for (j,q) in zip(eachindex(v.good_state),v.indices)
            (s[end-i+1,q]==v.good_state[end-i+1,j]) || return s, :undetected_failure
        end
    end
    return s, :true_success
end

"""Run a single Monte Carlo sample, starting with (and modifying) `initialstate` by applying the given `circuit`. Uses `applyop!` under the hood."""
function mctrajectory!(initialstate,circuit)
    state = initialstate
    for op in circuit
        state, cont = applyop!(state, op)
        if cont!=:continue
            return state, cont
        end
    end
    return state, :continue
end

"""Run multiple Monte Carlo trajectories and report the aggregate final statuses of each."""
function mctrajectories(initialstate,circuit;trajectories=500)
    counts = countmap([mctrajectory!(copy(initialstate),circuit)[2] for i in 1:trajectories]) # TODO use threads or at least a generator
    return merge(Dict([(k=>0) for k in statuses[2:end]]), counts)
end

"""Compute all possible new states after the application of the given operator. Reports the probability of each one of them. Deterministic, part of the Perturbative Expansion interface."""
function applyop_branches end

applyop_branches(s::Stabilizer, g::SparseGate; max_order=1) = [(applyop!(copy(s),g)...,1,0)] # there are no fall backs on purpose, otherwise it is easy to mistakenly make a non-deterministic version of this method
applyop_branches(s::Stabilizer, g::DenseGate; max_order=1) = [(applyop!(copy(s),g)...,1,0)]
applyop_branches(s::Stabilizer, v::VerifyOp; max_order=1) = [(applyop!(copy(s),v)...,1,0)] 

"""Compute all possible new states after the application of the given noise model. Reports the probability of each one of them. Deterministic, part of the Noise interface."""
function applynoise_branches end

function applynoise_branches(s::Stabilizer,noise::UnbiasedUncorrelatedNoise,indices::AbstractVector{Int}; max_order=1)
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
        for error_indices in combinations(indices, order) # TODO clean this up, optimize it
            for error_types in Base.Iterators.product(repeat([[apply_single_x!,apply_single_y!,apply_single_z!]],order)...)
                new_state = copy(s)
                for (i,t) in zip(error_indices, error_types)
                    t(new_state,i)
                end
                push!(results,(new_state, error_prob, order))
            end
        end
    end
    results
end

function applyop_branches(s::Stabilizer, nop::NoiseOpAll; max_order=1)
    n = nqubits(s)
    return [(state, :continue, prob, order) for (state, prob, order) in applynoise_branches(s, nop.noise, 1:n, max_order=max_order)]
end

function applyop_branches(s::Stabilizer, nop::NoiseOp; max_order=1)
    return [(state, :continue, prob, order) for (state, prob, order) in applynoise_branches(s, nop.noise, affectedqubits(nop), max_order=max_order)]
end

function applyop_branches(s::Stabilizer, g::NoisyGate; max_order=1)
    news, _,_,_ = applyop_branches(s,g.gate,max_order=max_order)[1] # TODO this assumes only one always successful branch for the gate
    return [(state, :continue, prob, order) for (state, prob, order) in applynoise_branches(news, g.noise, affectedqubits(g), max_order=max_order)]
end

function applyop_branches(s::Stabilizer, m::NoisyBellMeasurement; max_order=1)
    measurement_branches = applyop_branches(s, m.meas, max_order=max_order)
    if max_order==0
        return measurement_branches
    else
        new_branches = []
        nqubits = length(affectedqubits(m))
        errprob = 1//2*(1-p)
        sucprob = 1//2*(1+p)
        for (mstate, success, mprob, morder) in measurement_branches
            push!(new_branches, (mstate, success, mprob*sucprob, morder))
            push!(new_branches, (mstate, success==:continue ? :detected_failure : :continue, mprob*errprob, morder+1))
        end
        return new_branches
    end
end

# TODO a lot of repetition with applyop!
function applyop_branches(s::Stabilizer, m::BellMeasurement; max_order=1)
    n = nqubits(s)
    [(ns,iseven(r>>1) ? :continue : :detected_failure, p,0)
     for (ns,r,p) in _applyop_branches_measurement([(s,0x0,1.0)],m.pauli,affectedqubits(m),n)]
end

# TODO XXX THIS IS PARTICULARLY INEFFICIENT recurrent implementation
function _applyop_branches_measurement(branches, paulis, qubits, n)
    if length(paulis) == 0
        return branches
    end

    new_branches = []
    pauli = paulis[1]
    otherpaulis = paulis[2:end]
    index = qubits[1]
    otherqubits = qubits[2:end]
    if pauli == X # TODO this is not an elegant way to choose between X and Z coincidence measurements
        op = single_x(n,index) # TODO this is pretty terribly inefficient... use some sparse check
    elseif pauli == Z
        op = single_z(n,index)
    elseif pauli == Y
        op = single_y(n,index)
    else
        op = -single_y(n, index)
    end # TODO permit Y operators and permit negative operators

    for (s,r0,p) in branches
        s,anticom,r = project!(s,op)
        if isnothing(r)
            s1 = s
            s2 = copy(s)
            r1 = s1.phases[anticom] = 0x00
            r2 = s2.phases[anticom] = 0x02
            push!(new_branches, (s1,r0+r1,p/2))
            push!(new_branches, (s2,r0+r2,p/2))
        else
            push!(new_branches, (s,r0+r,p))
        end
    end

    return _applyop_branches_measurement(new_branches, otherpaulis, otherqubits, n)
end

# TODO a lot of repetition with applyop!
function applyop_branches(s::Stabilizer, mr::BellMeasurementAndReset; max_order=1)
    branches = applyop_branches(s,mr.meas, max_order=max_order)
    s = branches[1][1] # relies on the order of the branches, does not reset the branch with success==false, assumes order=0
    branches = [(_reset!(s,affectedqubits(mr).mr.resetto),succ,prob,order) for (s,succ,prob,order) in branches]
    branches
end

function _reset!(s, qubits, resetto)
    # TODO is the traceout necessary given that we just performed measurements?
    traceout!(s,qubits)# TODO it seems like a bad idea not to keep track of the rank here
    for (ii,i) in enumerate(qubits)
        for j in [1,2]
            s[end-j+1,i] = resetto[j,ii]
        end
    end
    return s
end

"""Run a perturbative expansion to a given order. Uses applyop_branches under the hood."""
function petrajectory(state, circuit; branch_weight=1.0, current_order=0, max_order=1)
    if size(circuit)[1] == 0
        return fill(zero(branch_weight), length(statuses)-1)
    end
    next_op = circuit[1]
    rest_of_circuit = circuit[2:end]

    status_probs = fill(zero(branch_weight), length(statuses)-1)

    for (i,(newstate, status, prob, order)) in enumerate(applyop_branches(state, next_op, max_order=max_order-current_order))
        if status==:continue # TODO is the copy below necessary?
            out_probs = petrajectory(copy(newstate), rest_of_circuit,
                branch_weight=branch_weight*prob, current_order=current_order+order, max_order=max_order)
            status_probs .+= out_probs
        else
            status_probs[findfirst(==(status),statuses)-1] += prob*branch_weight # TODO this findfirst needs to go, this is a bad way to do it
        end
    end

    return status_probs
end

"""Run a perturbative expansion to a given order. This is the main public fuction for the perturbative expansion approach."""
function petrajectories(state, circuit; branch_weight=1.0, max_order=1)
    status_probs = petrajectory(state, circuit; branch_weight=branch_weight, current_order=0, max_order=max_order)
    Dict([statuses[i+1]=>status_probs[i] for i in eachindex(status_probs)])
end

"""A register, representing the state of a computer including both a tableaux and an array of classical bits (e.g. for storing measurement results)"""
struct Register{S<:AbstractStabilizer, T<:AbstractVector{Bool}}
    stab::S
    bits::T
end

Base.copy(state::Register) = deepcopy(state)

function applyop!(state::Register, op)
    s, status = applyop!(state.stab, op)
    state, status
end

function applynoise!(state::Register, noise, indices)
    s, status = applynoise!(state.stab, noise, indices)
    state, status
end

function applyop_branches(state::Register, op; max_order=1)
    [(Register(newstate,copy(state.bits)), status, prob, order)
     for (newstate, status, prob, order)
     in applyop_branches(state.stab, op; max_order=max_order)
    ]
end

function applynoise_branches(state::Register, noise, indices; max_order=1)
    [(Register(newstate,copy(state.bits)), prob, order)
     for (newstate, prob, order) in applynoise_branches(s, nop.noise, 1:n, max_order=max_order)]
end

function applyop!(state::Register, op::DenseMeasurement)
    stab = state.stab
    stab,anticom,r = project!(stab, op.pauli)
    if isnothing(r)
        if rand()>0.5 # TODO float not necessary
            r = stab.phases[anticom] = 0x00
        else
            r = stab.phases[anticom] = 0x02
        end
    end
    state.bits[op.storagebit] = r==0x02
    state, :continue
end

function applyop!(state::Register, op::ConditionalGate)
    if state.bits[op.controlbit]
        applyop!(state, op.truegate)
    else
        applyop!(state, op.falsegate)
    end
    return state, :continue
end

function applyop!(state::Register, op::DecisionGate)
    decision = op.decisionfunction(state.bits)
    if !isnothing(decision)
        applyop!(state, op.gates[decision])
    end
    state, :continue
end


applyop_branches(s::Register, op::ConditionalGate; max_order=1) = [(applyop!(copy(s),op)...,1,0)]
applyop_branches(s::Register, op::DecisionGate; max_order=1) = [(applyop!(copy(s),op)...,1,0)]

function applyop_branches(s::Register, op::DenseMeasurement; max_order=1)
    stab = s.stab
    stab, anticom, r = project!(stab, op.pauli)
    new_branches = []
    if isnothing(r)
        s1 = s
        s1.stab.phases[anticom] = 0x00
        s1.bits[op.storagebit] = false
        s2 = copy(s)
        s2.stab.phases[anticom] = 0x02
        s2.bits[op.storagebit] = true
        push!(new_branches, (s1,:continue,1/2,0))
        push!(new_branches, (s2,:continue,1/2,0))
    else
        s.bits[op.storagebit] = r==0x02
        push!(new_branches, (s,:continue,1,0))
    end
    new_branches
end

include("./quantikz_methods.jl")

end
