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

end
