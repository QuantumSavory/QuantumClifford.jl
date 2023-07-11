#=
1. adding tests for basic correctness
2. single qubit gates / channels (and tests)
3. embedding single qubit gates - Stefan
4. pretty printing - Stefan
5. good docstrings
6. some superficial documentation
7. picking names
8. conversion into density matrices (QuantumOptics.jl) - Stefan
9. special small gates
10. make an overleaf for a paper
=#

"""
$(TYPEDEF)

Represents mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is a pure stabilizer state.

```jldoctest
julia> StabMixture(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> tT
Pauli channel ρ ↦ ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† with the following branches:
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.853553+0.0im | + _ | + _
 0.0-0.353553im | + _ | + Z
 0.0+0.353553im | + Z | + _
 0.146447+0.0im | + Z | + Z

julia> apply!(StabMixture(S"-X"), tT)
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.0-0.353553im | + _ | + Z
 0.0+0.353553im | + Z | + _
 0.853553+0.0im | + _ | + _
 0.146447+0.0im | + Z | + Z
```

See also: [`PauliChannel`](@ref)
"""
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
    print(io, "with ϕᵢⱼ | Pᵢ | Pⱼ:")
    for ((di,dj), χ) in s.destabweights
        println(io)
        print(io, " ")
        print(IOContext(io, :compact => true), χ)
        print(io, " | ", string(_stabmixdestab(s.stab, di)), " | ", string(_stabmixdestab(s.stab, dj)))
    end
end

function _stabmixdestab(mixeddestab, d)
    destab = destabilizerview(mixeddestab)
    p = zero(PauliOperator, nqubits(mixeddestab))
    for i in eachindex(d)
        if d[i]
            mul_right!(p, destab, i) # TODO check whether you are screwing up the ordering
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

"""A Pauli channel datastructure, mainly for use with [`StabMixture`](@ref)"""
struct PauliChannel{T,S} <: AbstractPauliChannel
    paulis::T
    weights::S
    function PauliChannel(paulis, weights)
        n = nqubits(paulis[1][1])
        for p in paulis
            n == nqubits(p[1]) == nqubits(p[2]) || throw(ArgumentError(lazy"""
            You are attempting to construct a `PauliChannel` but have provided Pauli operators
            $(p[1]) and $(p[2])
            that are not all of the same size (same number of qubits).
            Please ensure that all of the Pauli operators being provided of of the same size.
            """))
        end
        new{typeof(paulis),typeof(weights)}(paulis,weights)
    end
end

function Base.show(io::IO, pc::PauliChannel)
    println(io, "Pauli channel ρ ↦ ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† with the following branches:")
    print(io, "with ϕᵢⱼ | Pᵢ | Pⱼ:")
    for ((di,dj), χ) in zip(pc.paulis, pc.weights)
        println(io)
        print(io, " ")
        print(IOContext(io, :compact => true), χ)
        print(io, " | ", di, " | ", dj)
    end
end

function embed(n,idx,pc::PauliChannel)
    PauliChannel(map(p->(embed(n,idx,p[1]),embed(n,idx,p[2])),pc.paulis), pc.weights)
end

nqubits(pc::PauliChannel) = nqubits(pc.paulis[1][1])

function apply!(state::StabMixture, gate::PauliChannel)
    dict = state.destabweights
    stab = state.stab
    tzero = zero(eltype(dict).parameters[2])
    tone = one(eltype(dict).parameters[2])
    newdict = typeof(dict)(tzero) # TODO jeez, this is ugly
    for ((dᵢ,dⱼ), χ) in dict # the state
        for ((Pₗ,Pᵣ), w) in zip(gate.paulis,gate.weights) # the channel
            phaseₗ, dₗ, dₗˢᵗᵃᵇ = rowdecompose(Pₗ,stab)
            phaseᵣ, dᵣ, dᵣˢᵗᵃᵇ = rowdecompose(Pᵣ,stab)
            c = (dot(dₗˢᵗᵃᵇ,dᵢ) + dot(dᵣˢᵗᵃᵇ,dⱼ))*2
            dᵢ′ = dₗ .⊻ dᵢ
            dⱼ′ = dᵣ .⊻ dⱼ
            χ′ = χ * w * (-tone)^c * (tone*im)^(-phaseₗ+phaseᵣ+4)
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

"""Decompose a Pauli ``P`` in terms of stabilizer and destabilizer rows from a given tableaux.

For given tableaux of rows destabilizer rows ``\\{d_i\\}`` and stabilizer rows ``\\{s_i\\}``,
there are boolean vectors ``b`` and ``c`` such that
``P = i^p \\prod_i d_i^{b_i} \\prod_i s_i^{c_i}``.

This function returns `p`, `b`, `c`.

```
julia> s = MixedDestabilizer(ghz(2))
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _X
𝒮𝓉𝒶𝒷
+ XX
+ ZZ

julia> phase, destab_rows, stab_rows = QuantumClifford.rowdecompose(P"XY", s)
(3, Bool[1, 0], Bool[1, 1])

julia> im^3 * P"Z_" * P"XX" * P"ZZ"
+ XY
```
"""
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
            mul_right!(Pₜ, dest, i)
        end
    end
    for i in 1:n
        if comm(pauli, dest, i) != 0
            c[i] = true
            mul_right!(Pₜ, stab, i)
        end
    end
    p = mod(-Pₜ.phase[],4) # complex conjugate
    return p, b, c
end

##
# Predefined Pauli Channels
##

const tT = PauliChannel(
    [(I, I), (I, Z), (Z, I), (Z, Z)],
    [cos(π/8)^2, -im*sin(π/8)*cos(π/8),  im*sin(π/8)*cos(π/8), sin(π/8)^2]
    )
