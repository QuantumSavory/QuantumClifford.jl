module NoisyCircuits

# TODO the current interfaces of this module are a bit clunky when it comes to random events
# - applyop! returns (new_state, confirmation_of_no_detected_errors)
# - two sources of randomness are not tracked well currently
#   - which noise operator is applied (applied in a purely MC fashion)
#   - when a measurement is not commuting, which branch is taken (purely MC for now as well)
#
# TODO should we call things `apply!` instead of applyop and applynoise?
#
# TODO how important it is to distinguish measuring X₁ and then X₂ from measuring X₁X₂ when doing coincidence measurements

using QuantumClifford

export Operation, AbstractGate, AbstractMeasurement, AbstractNoise,
       UnbiasedUncorrelatedNoise, NoiseOp, NoiseOpAll,
       SparseGate, NoisyGate,
       Measurement, NoisyMeasurement, MeasurementAndReset, MeasurementAndNoisyReset,
       affectedqubits, applyop!, applynoise!,
       NoisyCircuitResult, undetected_failure, detected_failure, true_success, mctrajectory!,
       petrajectory

abstract type Operation end
abstract type AbstractGate <: Operation end
abstract type AbstractMeasurement <: Operation end

abstract type AbstractNoise end

struct UnbiasedUncorrelatedNoise <: AbstractNoise
    errprobthird::Float64 # TODO should a float be hardcoded?
end

struct NoiseOp <: Operation
    noise::AbstractNoise
    indices::AbstractVector{Int}
end

struct NoiseOpAll <: Operation
    noise::AbstractNoise
end

struct SparseGate <: AbstractGate
    cliff::CliffordOperator # TODO do not hardcode this type of clifford op
    indices::AbstractVector{Int}
end

struct NoisyGate <: AbstractGate
    gate::AbstractGate # TODO should the type be more specific
    noise::AbstractNoise # TODO should the type be more specific
end

struct Measurement <: AbstractMeasurement
    pauli::AbstractVector{PauliOperator}
    indices::AbstractVector{Int}
end

struct NoisyMeasurement <: AbstractMeasurement
    meas::AbstractMeasurement # TODO should the type be more specific
    noise::AbstractNoise # TODO should the type be more specific
end

struct MeasurementAndReset <: AbstractMeasurement # TODO Do we need a new type or should all non-terminal measurements implicitly have a reset?
    meas::AbstractMeasurement # TODO is this the cleanest way to specify the type
    resetto::Stabilizer
end

struct MeasurementAndNoisyReset <: AbstractMeasurement # TODO Do we need a new type or should all non-terminal measurements implicitly have a reset?
    meas::AbstractMeasurement # TODO is this the cleanest way to specify the type
    resetto::Stabilizer
    noise::AbstractNoise # TODO should the type be more specific
end

affectedqubits(g::NoisyGate) = affectedqubits(g.gate)
affectedqubits(g::SparseGate) = g.indices
affectedqubits(m::Measurement) = m.indices
affectedqubits(m::MeasurementAndReset) = affectedqubits(m.meas)
affectedqubits(m::MeasurementAndNoisyReset) = affectedqubits(m.meas)
affectedqubits(m::NoisyMeasurement) = affectedqubits(m.meas)
affectedqubits(m::NoiseOp) = m.indices

function applyop!(s::Stabilizer, g::NoisyGate)
    s = applynoise!(
            applyop!(s,g.gate)[1],
            g.noise,
            affectedqubits(g.gate)),
    return s, true
end

applyop!(s::Stabilizer, g::SparseGate) = (apply!(s,g.cliff,affectedqubits(g)), true)

applyop!(s::Stabilizer, m::NoisyMeasurement) =  applyop!(
    applynoise!(s,m.noise,affectedqubits(m)),
    m.meas)

# TODO this seems unnecessarily complicated
function applyop!(s::Stabilizer, m::Measurement) # TODO is it ok to just measure XX instead of measuring XI and IX separately? That would be much faster
    n = nqubits(s)
    indices = affectedqubits(m)
    res = 0x00
    for (pauli, index) in zip(m.pauli,affectedqubits(m))
        if pauli==X # TODO this is not an elegant way to choose between X and Z coincidence measurements
            op = single_x(n,index) # TODO this is pretty terribly inefficient... use some sparse check
        else
            op = single_z(n,index)
        end # TODO permit Y operators and permit negative operators
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
        return s, true
    else
        return s, false
    end
end

function applyop!(s::Stabilizer, mr::MeasurementAndReset)
    s,res = applyop!(s,mr.meas)
    if !res
        return s,res
    else
        # TODO is the traceout necessary given that we just performed measurements?
        traceout!(s,mr.meas.indices)# TODO it seems like a bad idea not to keep track of the rank here
        n = nqubits(s) # TODO implement lastindex so we can just use end
        for (ii,i) in enumerate(affectedqubits(mr))
            for j in [1,2]
                s[n-j+1,i] = mr.resetto[j,ii]
            end
        end
        return s,true
    end
end

function applyop!(s::Stabilizer, mr::MeasurementAndNoisyReset)
    s,res = applyop!(s,mr.meas)
    if !res
        return s,res
    else
        # TODO is the traceout necessary given that we just performed measurements?
        traceout!(s,affectedqubits(mr))# TODO it seems like a bad idea not to keep track of the rank here
        n = nqubits(s) # TODO implement lastindex so we can just use end
        for (ii,i) in enumerate(affectedqubits(mr))
            for j in [1,2]
                s[n-j+1,i] = mr.resetto[j,ii]
            end
        end
        return applynoise!(s,mr.noise,affectedqubits(mr)), true
    end
end

function applynoise!(s::Stabilizer,noise::UnbiasedUncorrelatedNoise,indices::AbstractVector{Int})
    n = nqubits(s)
    infid = noise.errprobthird
    for i in indices
        r = rand()
        if r<infid
            apply!(s,single_x(n,i)) # TODO stupidly inefficient, do it sparsely
        end
        if infid<=r<2infid
            apply!(s,single_z(n,i)) # TODO stupidly inefficient, do it sparsely
        end
        if 2infid<=r<3infid
            apply!(s,single_x(n,i)) # TODO stupidly inefficient, do it sparsely
            apply!(s,single_z(n,i)) # TODO stupidly inefficient, do it sparsely
        end
    end
    s
end

function applyop!(s::Stabilizer, mr::NoiseOpAll)
    n = nqubits(s)
    return applynoise!(s, mr.noise, 1:n), true
end

function applyop!(s::Stabilizer, mr::NoiseOp)
    return applynoise!(s, mr.noise, affectedqubits(mr)), true
end

@enum NoisyCircuitResult undetected_failure detected_failure true_success

function mctrajectory!(initialstate::Stabilizer,circuit::AbstractVector{Operation},is_good)
    state = initialstate
    for op in circuit
        #println(typeof(op))
        state, success = applyop!(state, op)
        #println("#",typeof(state))
        if !success
            return state, detected_failure
        end
    end
    if is_good(state)
        return state, true_success
    else
        return state, undetected_failure
    end
end

applyop_branches(s::Stabilizer, g::SparseGate; max_order=1) = [(applyop!(copy(s),g)...,1,0)]

function applynoise_branches(s::Stabilizer,noise::UnbiasedUncorrelatedNoise,indices::AbstractVector{Int}; max_order=1)
    n = nqubits(s)
    l = length(indices)
    infid = noise.errprobthird
    if l==0
        return [s,one(infid)]
    end
    no_error1 = 1-3*infid
    no_error = no_error1^l
    single_error = no_error1^(l-1)*infid
    results = [(copy(s),no_error,0)] # state, prob, order
    if max_order==0
        return results
    end
    for i in indices # TODO max_order>1 is not currently implemented
        push!(results,(apply!(copy(s),single_x(n,i)), single_error, 1)) # TODO stupidly inefficient, do it sparsely
        push!(results,(apply!(copy(s),single_z(n,i)), single_error, 1)) # TODO stupidly inefficient, do it sparsely
        push!(results,(apply!(apply!(copy(s),single_x(n,i)),single_z(n,i)), single_error, 1)) # TODO stupidly inefficient, do it sparsely
    end
    results
end

function applyop_branches(s::Stabilizer, nop::NoiseOpAll; max_order=1)
    n = nqubits(s)
    return [(state, true, prob, order) for (state, prob, order) in applynoise_branches(s, nop.noise, 1:n, max_order=max_order)]
end

function applyop_branches(s::Stabilizer, nop::NoiseOp; max_order=1)
    return [(state, true, prob, order) for (state, prob, order) in applynoise_branches(s, nop.noise, affectedqubits(nop), max_order=max_order)]
end

function applyop_branches(s::Stabilizer, g::NoisyGate; max_order=1)
    news, _,_,_ = applyop_branches(s,g.gate,max_order=max_order)[1] # TODO this assumes only one always successful branch for the gate
    return [(state, true, prob, order) for (state, prob, order) in applynoise_branches(news, g.noise, affectedqubits(g), max_order=max_order)]
end

# TODO this can be much faster if we perform the flip on the classical bit after measurement, when possible
function applyop_branches(s::Stabilizer, m::NoisyMeasurement; max_order=1)
    return [(state, success, nprob*mprob, order)
            for (mstate, success, mprob, morder) in applyop_branches(s, m.meas, max_order=max_order)
            for (state, nprob, order) in applynoise_branches(mstate, m.noise, affectedqubits(m), max_order=max_order-morder)]
end

# TODO a lot of repetition with applyop!
function applyop_branches(s::Stabilizer, m::Measurement; max_order=1) # TODO is it ok to just measure XX instead of measuring XI and IX separately? That would be much faster
    n = nqubits(s)
    [(ns,iseven(r>>1),p,0) for (ns,r,p) in _applyop_branches_measurement([(s,0x0,1.0)],m.pauli,affectedqubits(m),n)]
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
    if pauli==X # TODO this is not an elegant way to choose between X and Z coincidence measurements
        op = single_x(n,index) # TODO this is pretty terribly inefficient... use some sparse check
    else
        op = single_z(n,index)
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
function applyop_branches(s::Stabilizer, mr::MeasurementAndReset; max_order=1)
    branches = applyop_branches(s,mr.meas, max_order=max_order)
    s = branches[1][1] # relies on the order of the branches, does not reset the branch with success==false, assumes order=0
    branches = [(_reset!(s,affectedqubits(mr).mr.resetto),succ,prob,order) for (s,succ,prob,order) in branches]
    branches
end

# TODO a lot of repetition with applyop!
function applyop_branches(s::Stabilizer, mr::MeasurementAndNoisyReset; max_order=1)
    branches = applyop_branches(s,mr.meas)
    # TODO can skip the inner loop if succ=false
    noise_branches = [(state, succ, prob*mprob, order) 
                      for (ms, succ, mprob, prime_order) in branches
                      for (state, prob, order) in applynoise_branches(_reset!(ms,affectedqubits(mr),mr.resetto), mr.noise, affectedqubits(mr), max_order=max_order-prime_order)
                     ]
end

function _reset!(s, qubits, resetto)
    # TODO is the traceout necessary given that we just performed measurements?
    traceout!(s,qubits)# TODO it seems like a bad idea not to keep track of the rank here
    n = nqubits(s) # TODO implement lastindex so we can just use end
    for (ii,i) in enumerate(qubits)
        for j in [1,2]
            s[n-j+1,i] = resetto[j,ii]
        end
    end
    return s
end

function petrajectory(state, circuit, is_good; branch_weight=1.0, current_order=0, max_order=1)
    if length(circuit)==0 # end case of the recursion
        if is_good(state)
            return (undetected_failure=zero(branch_weight),
            detected_failure=zero(branch_weight),
            true_success=branch_weight)
        else
            return (undetected_failure=branch_weight,
                    detected_failure=zero(branch_weight),
                    true_success=zero(branch_weight))
        end
    end

    next_op = circuit[1]
    rest_of_circuit = circuit[2:end]

    P_undetected_failure, P_detected_failure, P_true_success = zero(branch_weight), zero(branch_weight), zero(branch_weight)

    # applyop_all returns all branches of the noise model
    p = 0
    for (i,(newstate, success, prob, order)) in enumerate(applyop_branches(state, next_op, max_order=max_order-current_order))
        p+=prob
        if success # TODO is the copy below necessary?
            uf,df,ts = petrajectory(copy(newstate), rest_of_circuit, is_good,
                branch_weight=branch_weight*prob, current_order=current_order+order, max_order=max_order)
            P_undetected_failure += uf
            P_detected_failure += df
            P_true_success += ts
        else
            P_detected_failure += prob*branch_weight
        end
    end

    return (undetected_failure=P_undetected_failure,
    detected_failure=P_detected_failure,
    true_success=P_true_success)
end

end
