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
       NoisyCircuitResult, undetected_failure, detected_failure, true_success, mctrajectory!

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
        return applynoise!(s,mr.noise,affectedqubits(mr)), true # TODO indices should be a multimethod, not a property
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
    for i in indices
        push!((apply!(copy(s),single_x(n,i)), single_error, 1)) # TODO stupidly inefficient, do it sparsely
        push!((apply!(copy(s),single_z(n,i)), single_error, 1)) # TODO stupidly inefficient, do it sparsely
        push!((apply!(apply!(copy(s),single_x(n,i)),single_z(n,i)), single_error, 1)) # TODO stupidly inefficient, do it sparsely
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
    return [(state, true, prob, order) for (state, prob, order) in applynoise_branches(s, g.noise, affectedqubits(g), max_order=max_order)]
end

function applyop_branches(s::Stabilizer, m::NoisyMeasurement; max_order=1)
    return [(state, success, nprob*mprob, order)
            for (mstate, success, mprob, morder) in applyop_branches(s, g.meas, max_order=max_order)
            for (state, nprob, order) in applynoise_branches(mstate, g.noise, affectedqubits(g), max_order=max_order-morder)]
end

# TODO a lot of repetition with applyop!
# TODO XXX THIS IS COMPLETELY UNFINISHED - it performs only the first necessary measurement
function applyop_branches(s::Stabilizer, m::Measurement; max_order=1) # TODO is it ok to just measure XX instead of measuring XI and IX separately? That would be much faster
    n = nqubits(s)
    indices = affectedqubits(m)
    for (pauli, index) in zip(m.pauli,affectedqubits(m))
        if pauli==X # TODO this is not an elegant way to choose between X and Z coincidence measurements
            op = single_x(n,index) # TODO this is pretty terribly inefficient... use some sparse check
        else
            op = single_z(n,index)
        end # TODO permit Y operators and permit negative operators
        s,anticom,r = project!(copy(s),op)
        if isnothing(r)
            s1 = s
            s2 = copy(s)
            r1 = s1.phases[anticom] = 0x00
            r2 = s2.phases[anticom] = 0x02
            return [(s1,true,0.5,0), (s2,false,0.5,0)]
        else
            return [(s,r==0x0,1,0)]
        end
    end
end

# TODO a lot of repetition with applyop!
function applyop_branches(s::Stabilizer, mr::MeasurementAndReset; max_order=1)
    branches = applyop_branches(s,mr.meas, max_order=max_order)
    s = branches[1][1] # relies on the order of the branches, does not reset the branch with success==false, assumes order=0
    # TODO is the traceout necessary given that we just performed measurements?
    traceout!(s,mr.meas.indices)# TODO it seems like a bad idea not to keep track of the rank here
    n = nqubits(s) # TODO implement lastindex so we can just use end
    for (ii,i) in enumerate(affectedqubits(mr))
        for j in [1,2]
            s[n-j+1,i] = mr.resetto[j,ii]
        end
    end
    return branches
end

# TODO a lot of repetition with applyop!
function applyop_branches(s::Stabilizer, mr::MeasurementAndNoisyReset; max_order=1)
    branches = applyop_branches(s,mr.meas)
    ms, _, mprob, _ =  branches[1]# relies on the order of the branches, does not reset the branch with success==false, assumes order=0
    noise_branches = [(state, true, prob*mprob, order) for (state, prob, order) in applynoise_branches(ms, mr.noise, affectedqubits(ms), max_order=max_order)]
    return vcat(noise_branches,branches[2:end])
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
    for (newstate, success, prob, order) in applyop_branches(state, next_op, max_order=max_order-current_order)
        println("-"^order," ", prob)
        if success
            uf,df,ts = petrajectory(newstate, rest_of_circuit, is_good,
                branch_weight=branch_weight*prob, current_order=current_order+order, max_order=max_order)
            P_undetected_failure += uf
            P_detected_failure += df
            P_true_success += ts
        else
            P_detected_failure += prob
        end
    end

    return (undetected_failure=P_undetected_failure,
            detected_failure=P_detected_failure,
            true_success=P_true_success)
end

end
