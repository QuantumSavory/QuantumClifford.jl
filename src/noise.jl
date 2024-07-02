abstract type AbstractNoise end

abstract type AbstractNoiseOp <: AbstractOperation end

"""A method modifying a given state by applying the corresponding noise model. It is non-deterministic, part of the Noise interface."""
function applynoise! end

function applynoise!(state, noise, indices::Base.AbstractVecOrTuple)
    @inbounds @simd for i in indices
        @inline applynoise!(state, noise, i)
    end
    return state
end

# Implementations for Register
function applynoise!(r::Register, n, i::Int)
    apply!(quantumstate(r), n, i)
    return r
end
function applynoise!(r::Register, n, indices::Base.AbstractVecOrTuple)
    apply!(quantumstate(r), n, indices)
    return r
end

"""Depolarization noise model with total probability of error `p`."""
struct UnbiasedUncorrelatedNoise{T} <: AbstractNoise
    p::T
end
UnbiasedUncorrelatedNoise(p::Integer) = UnbiasedUncorrelatedNoise(float(p))

"""Pauli noise model with probabilities `px`, `py`, and `pz` respectively for the three types of Pauli errors."""
struct PauliNoise{T} <: AbstractNoise
    px::T
    py::T
    pz::T
end
function PauliNoise(px::Real, py::Real, pz::Real)
    px, py, pz = float.((px, py, pz))
    px, py, pz = promote(px, py, pz)
    T = typeof(px)
    return PauliNoise{T}(px, py, pz)
end

"""A convenient constructor for various types of Pauli noise models.
Returns more specific types when necessary."""
function PauliNoise end

"""Constructs an unbiased Pauli noise model with total probability of error `p`."""
function PauliNoise(p)
    UnbiasedUncorrelatedNoise(p)
end

function applynoise!(s::AbstractStabilizer,noise::UnbiasedUncorrelatedNoise,i::Int)
    infid = noise.p/3
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

function applynoise!(s::AbstractStabilizer,noise::PauliNoise,i::Int)
    r = rand()
    if r<noise.px
        apply_single_x!(s,i)
    elseif r<noise.px+noise.pz
        apply_single_z!(s,i)
    elseif r<noise.px+noise.pz+noise.py
        apply_single_y!(s,i)
    end
    s
end

"""An operator that applies the given `noise` model to the qubits at the selected `indices`."""
struct NoiseOp{N, Q} <: AbstractNoiseOp where {N, Q}
    noise::N #<:AbstractNoise
    indices::NTuple{Q, Int}
end

NoiseOp(noise, indices::AbstractVector{Int}) = NoiseOp(noise, tuple(indices...))

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

function PauliError(q::Int, px, py, pz)
    NoiseOp(PauliNoise(px,py,pz), (q,))
end

function PauliError(qubits, px, py, pz)
    NoiseOp(PauliNoise(px,py,pz), qubits)
end

"""An operator that applies the given `noise` model to all qubits."""
struct NoiseOpAll <: AbstractNoiseOp
    noise::AbstractNoise
end

"""A gate consisting of the given `noise` applied after the given perfect Clifford `gate`."""
struct NoisyGate <: AbstractNoiseOp
    gate::AbstractOperation
    noise::AbstractNoise
end

function apply!(s::AbstractQCState, g::NoisyGate)
    s = applynoise!(
            apply!(s,g.gate),
            g.noise,
            affectedqubits(g.gate)),
    return s
end

function apply!(s::AbstractQCState, mr::NoiseOpAll)
    n = nqubits(s)
    return applynoise!(s, mr.noise, 1:n)
end

function apply!(s::AbstractQCState, mr::NoiseOp)
    return applynoise!(s, mr.noise, affectedqubits(mr))
end

# XXX necessary to resolve ambiguity between apply!(s::AbstractQCState, mr::Noise) and apply!(r::Register, op)
# TODO resolve them in a neater fashion with less repetition
function apply!(r::Register, n::NoisyGate)
    apply!(quantumstate(r), n)
    return r
end
function apply!(r::Register, n::NoiseOpAll)
    apply!(quantumstate(r), n)
    return r
end
function apply!(r::Register, n::NoiseOp)
    apply!(quantumstate(r), n)
    return r
end

##
# petrajectories
##

function applynoise_branches(s::AbstractStabilizer,noise::UnbiasedUncorrelatedNoise,indices; max_order=1)
    n = nqubits(s)
    l = length(indices)
    infid = noise.p/3
    if l==0
        return [s,one(infid)]
    end
    error1 = noise.p
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
