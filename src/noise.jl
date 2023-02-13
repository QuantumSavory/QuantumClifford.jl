abstract type AbstractNoise end

"""A method modifying a given state by applying the corresponding noise model. Non-deterministic, part of the Noise interface."""
function applynoise! end

function applynoise!(state::Register, noise, indices)
    s = applynoise!(state.stab, noise, indices)
    state
end

"""Depolarization noise model with total probability of error `3*errprobthird`."""
struct UnbiasedUncorrelatedNoise{T} <: AbstractNoise
    errprobthird::T
end

function applynoise!(s::AbstractStabilizer,noise::UnbiasedUncorrelatedNoise,indices::AbstractVector{Int})
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

"""An operator that applies the given `noise` model to the qubits at the selected `indices`."""
struct NoiseOp <: AbstractOperation
    noise::AbstractNoise
    indices::AbstractVector{Int}
end

"""An operator that applies the given `noise` model to all qubits."""
struct NoiseOpAll <: AbstractOperation
    noise::AbstractNoise
end

"""A gate consisting of the given `noise` applied after the given perfect Clifford `gate`."""
struct NoisyGate <: AbstractOperation
    gate::AbstractOperation
    noise::AbstractNoise
end

function apply!(s::AbstractStabilizer, g::NoisyGate)
    s = applynoise!(
            apply!(s,g.gate),
            g.noise,
            affectedqubits(g.gate)),
    return s
end

function apply!(s::AbstractStabilizer, mr::NoiseOpAll)
    n = nqubits(s)
    return applynoise!(s, mr.noise, 1:n)
end

function apply!(s::AbstractStabilizer, mr::NoiseOp)
    return applynoise!(s, mr.noise, affectedqubits(mr))
end
