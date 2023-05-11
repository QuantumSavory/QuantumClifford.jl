abstract type AbstractNoise end

"""A method modifying a given state by applying the corresponding noise model. It is non-deterministic, part of the Noise interface."""
function applynoise! end

function applynoise!(state, noise, indices)
    @inbounds @simd for i in indices
        @inline applynoise!(state, noise, i)
    end
end

"""Depolarization noise model with total probability of error `3*errprobthird`."""
struct UnbiasedUncorrelatedNoise{T} <: AbstractNoise
    errprobthird::T
end

"""A convenient constructor for various types of Pauli noise models.
Returns more specific types when necessary."""
function PauliNoise end

"""Constructs an unbiased Pauli noise model with total probability of error `p`."""
function PauliNoise(p)
    UnbiasedUncorrelatedNoise(p/3)
end

function applynoise!(s::AbstractStabilizer,noise::UnbiasedUncorrelatedNoise,i::Int)
    n = nqubits(s)
    infid = noise.errprobthird
    r = rand()
    if r<infid
        apply_single_x!(s,i)
    elseif r<2infid
        apply_single_z!(s,i)
    elseif r<3infid
        apply_single_y!(s,i)
    end
    s
end

"""An operator that applies the given `noise` model to the qubits at the selected `indices`."""
struct NoiseOp <: AbstractOperation
    noise::AbstractNoise
    indices::AbstractVector{Int}
end

"""A convenient constructor for various types of Pauli errors,
that can be used as circuit gates in simulations.
Returns more specific types when necessary."""
function PauliError end

""""Construct a gate operation that applies an unbiased Pauli error on qubit `q` with probability `p`."""
function PauliError(q::Int,p)
    NoiseOp(PauliNoise(p), [q])
end

""""Construct a gate operation that applies an unbiased Pauli error on all `qubits`, each with independent probability `p`."""
function PauliError(qubits,p)
    NoiseOp(PauliNoise(p), qubits)
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
