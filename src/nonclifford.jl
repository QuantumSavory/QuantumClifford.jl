using LinearAlgebra: dot
using DataStructures: DefaultDict

mutable struct StabMixture{T,F}
    stab::T
    destabweights::DefaultDict{Tuple{BitVector, BitVector}, F, F}
end

function StabMixture(state)
    n = nqubits(state)
    StabMixture(MixedDestabilizer(state), DefaultDict(0.0im, (falses(n),falses(n))=>1.0+0.0im)) # TODO maybe it should default to Destabilizer, not MixedDestabilizer
end

StabMixture(s::StabMixture) = s

function MixedDestabilizer(s::StabMixture)
    if length(s.destabweights) != 1
        throw(DomainError("Trying to convert a non-Clifford state (instance of type StabMixture with more than one term in its sum representation) to a Clifford state (of type MixedDestabilizer). This is not possible. Consider whether you have performed some (unforeseen) non-Clifford operation on the state."))
    else
        ((dᵢ, dⱼ), χ) = first(s.destabweights)
        dᵢ = _stabmixdestab(s.stab, dᵢ) # TODO
        dⱼ = _stabmixdestab(s.stab, dⱼ) # make
        stab = copy(s.stab)             # a faster
        mul_left!(stab, dᵢ)             # in-place
        mul_right!(stab, dⱼ)            # version
        return stab
    end
end

function Base.show(io::IO, s::StabMixture)
    println(io, "A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is")
    show(io,s.stab)
    println(io)
    println(io, "with ϕᵢⱼ | Pᵢ | Pⱼ:")
    for ((di,dj), χ) in s.destabweights
        println(io, "  ", χ, " | ", string(_stabmixdestab(s.stab, di)), " | ", string(_stabmixdestab(s.stab, dj)))
    end
end

function _stabmixdestab(mixeddestab, d)
    destab = destabilizerview(mixeddestab)
    p = zero(PauliOperator, nqubits(mixeddestab))
    for i in eachindex(d)
        if d[i]
            mul_left!(p, destab, i) # TODO check whether you are screwing up the ordering
        end
    end
    p
end

function apply!(state::StabMixture, gate::AbstractCliffordOperator) # TODO conjugate also the destabs
    apply!(state.stab, gate)
    state
end

abstract type AbstractPauliChannel <: AbstractOperation end

struct TGate <: AbstractPauliChannel
    qubit::Int
end

struct PauliChannel{T,S} <: AbstractPauliChannel
    paulis::T
    weights::S
end

function apply!(state::StabMixture, gate::PauliChannel)
    dict = state.destabweights
    stab = state.stab
    newdict = typeof(dict)(zero(eltype(dict).parameters[2])) # TODO jeez, this is ugly
    for ((dᵢ,dⱼ), χ) in dict
        for ((Pₗ,Pᵣ), w) in zip(gate.paulis,gate.weights)
            phaseₗ, dₗ, dₗˢᵗᵃᵇ = rowdecompose(Pₗ,stab)
            phaseᵣ, dᵣ, dᵣˢᵗᵃᵇ = rowdecompose(Pᵣ,stab)
            c = (dot(dₗˢᵗᵃᵇ,dᵢ) + dot(dᵣˢᵗᵃᵇ,dⱼ))*2
            dᵢ′ = dₗ .⊻ dᵢ
            dⱼ′ = dᵣ .⊻ dⱼ
            χ′ = χ * w * (-1)^c * (1im)^(phaseₗ-phaseᵣ)
            newdict[(dᵢ′,dⱼ′)] += χ′
        end
    end
    for (k,v) in newdict # TODO is it safe to modify a dict while iterating over it?
        if abs(v) < 1e-14 # TODO parameterize this pruning parameter
            delete!(newdict, k)
        end
    end
    state.destabweights = newdict
    state
end

function rowdecompose(pauli,state::Union{MixedDestabilizer, Destabilizer})
    n = nqubits(pauli)
    stab = stabilizerview(state)
    dest = destabilizerview(state)
    b = falses(n)
    c = falses(n)
    Pₜ = zero(PauliOperator, n) # TODO is this the best API choice!?
    for i in 1:n
        if comm(pauli, stab, i) != 0
            b[i] = true
            mul_left!(Pₜ, dest, i)
        end
    end
    for i in 1:n
        if comm(pauli, dest, i) != 0
            c[i] = true
            mul_left!(Pₜ, stab, i) # ?
        end
    end
    return Pₜ.phase[], b, c
end
