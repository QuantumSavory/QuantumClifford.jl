import QuantumClifford: ⊗

import Base: *, copy

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
julia> GeneralizedStabilizer(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> GeneralizedStabilizer(S"-X").destabweights
DataStructures.DefaultDict{Tuple{BitVector, BitVector}, ComplexF64, ComplexF64} with 1 entry:
  ([0], [0]) => 1.0+0.0im

julia> pcT
A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:
with ϕᵢ | Pᵢ
 0.853553+0.353553im | + _
 0.146447-0.353553im | + Z

julia> gs = apply!(GeneralizedStabilizer(S"-X"), pcT)
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.0+0.353553im | + _ | + Z
 0.0-0.353553im | + Z | + _
 0.853553+0.0im | + _ | + _
 0.146447+0.0im | + Z | + Z

julia> gs.stab
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X

julia> gs.destabweights
DataStructures.DefaultDict{Tuple{BitVector, BitVector}, ComplexF64, ComplexF64} with 4 entries:
  ([0], [1]) => 0.0+0.353553im
  ([1], [0]) => 0.0-0.353553im
  ([0], [0]) => 0.853553+0.0im
  ([1], [1]) => 0.146447+0.0im
```

See also: [`PauliChannel`](@ref)
"""
mutable struct GeneralizedStabilizer{T,F}
    stab::T
    destabweights::DefaultDict{Tuple{BitVector, BitVector}, F, F}
end

function GeneralizedStabilizer(state)
    n = nqubits(state)
    md = MixedDestabilizer(state)
    rank(md)==n || throw(ArgumentError(lazy"""
        Attempting to convert a `Stabilizer`-like object to `GeneralizedStabilizer` object failed,
        because the initial state does not represent a pure state.
        Currently only pure states can be used to initialize a `GeneralizedStabilizer` mixture of stabilizer states.
    """))
    GeneralizedStabilizer(md, DefaultDict(0.0im, (falses(n),falses(n))=>1.0+0.0im)) # TODO maybe it should default to Destabilizer, not MixedDestabilizer
end

GeneralizedStabilizer(s::GeneralizedStabilizer) = s

function Base.show(io::IO, s::GeneralizedStabilizer)
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

function apply!(state::GeneralizedStabilizer, gate::AbstractCliffordOperator) # TODO conjugate also the destabs
    apply!(state.stab, gate)
    state
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> sm = GeneralizedStabilizer(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> sm ⊗ sm
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _Z
𝒮𝓉𝒶𝒷
- X_
- _X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + __ | + __
```

"""
function (⊗)(state₁::GeneralizedStabilizer, state₂::GeneralizedStabilizer)
    dict₁ = state₁.destabweights
    dict₂ = state₂.destabweights
    dtype = valtype(dict₁)
    tzero = zero(dtype)
    newdict = typeof(dict₁)(tzero)
    newstab = state₁.stab ⊗ state₂.stab
    n = nqubits(newstab)
    newsm = GeneralizedStabilizer(newstab, DefaultDict(0.0im, (falses(n),falses(n))=>1.0+0.0im))
    for ((dᵢ, dⱼ), χ) in dict₁
        for ((dᵢ′, dⱼ′), χ′) in dict₂
            newdict[(dᵢ, dⱼ)] = get!(newdict, (dᵢ, dⱼ), tzero) + χ * χ′
        end
    end
    newsm.destabweights = newdict
    return newsm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> sm = GeneralizedStabilizer(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> tHadamard * sm
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ X
𝒮𝓉𝒶𝒷
- Z
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _
```

"""
function (*)(Op::AbstractCliffordOperator, state::GeneralizedStabilizer)
    dict = state.destabweights
    dtype = valtype(dict)
    tzero = zero(dtype)
    newdict = typeof(dict)(tzero)
    newstab = Op * state.stab
    n = nqubits(newstab)
    newsm = GeneralizedStabilizer(newstab, DefaultDict(0.0im, (falses(n),falses(n))=>1.0+0.0im))
    newsm.destabweights = dict
    return newsm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> sm = GeneralizedStabilizer(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> P"-Y" * sm
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
- Z
𝒮𝓉𝒶𝒷
+ X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _
```

"""
function (*)(Op::PauliOperator, state::GeneralizedStabilizer)
    dict = state.destabweights
    dtype = valtype(dict)
    tzero = zero(dtype)
    newdict = typeof(dict)(tzero)
    newstab = Op * state.stab
    n = nqubits(newstab)
    newsm = GeneralizedStabilizer(newstab, DefaultDict(0.0im, (falses(n),falses(n))=>1.0+0.0im))
    newsm.destabweights = dict
    return newsm
end 


"""$(TYPEDSIGNATURES)

Expectation value for the [PauliOperator](@ref) observable given the [`GeneralizedStabilizer`](@ref) state `s`."""
function expect(p::PauliOperator, s::GeneralizedStabilizer) # TODO optimize
    e = zero(_dictvaltype(s.destabweights))
    phase, b, c = rowdecompose(p, s.stab)
    for ((dᵢ,dⱼ), χ) in s.destabweights
        _allthreesumtozero(dᵢ,dⱼ,b) || continue
        e += χ * (-1)^(dᵢ'*c)
    end
    return (-1)^(phase÷2) * e
end

"""Same as `all(==(0), (a.+b.+c) .% 2)`"""
function _allthreesumtozero(a,b,c) # TODO consider using bitpacking and SIMD xor with less eager shortcircuiting -- probably would be much faster
    @inbounds @simd for i in 1:length(a)
        ((a[i] + b[i] + c[i]) & 1) == 0 || return false
    end
    true
end

function _dictvaltype(dict)
    return valtype(dict)
end

function project!(sm::GeneralizedStabilizer, p::PauliOperator)
    eval = expect(p, sm)
    prob₁ = (real(eval)+1)/2
    error("This functionality is not implemented yet")
end

function _proj₋(sm::GeneralizedStabilizer, p::PauliOperator)
end
function _proj₊(sm::GeneralizedStabilizer, p::PauliOperator)
end

abstract type AbstractPauliChannel <: AbstractOperation end

"""A Pauli channel datastructure, mainly for use with [`GeneralizedStabilizer`](@ref)

See also: [`UnitaryPauliChannel`](@ref)"""
struct PauliChannel{T,S} <: AbstractPauliChannel
    paulis::T
    weights::S
    function PauliChannel(paulis, weights)
        length(paulis) == length(weights) || throw(ArgumentError(lazy"""
        Attempting to construct a `PauliChannel` failed.
        The length of the vectors of weights and of Pauli operator pairs differs
        ($(length(weights)) and $(length(paulis)) respectively).
        """))
        n = nqubits(paulis[1][1])
        for p in paulis
            n == nqubits(p[1]) == nqubits(p[2]) || throw(ArgumentError(lazy"""
            You are attempting to construct a `PauliChannel` but have provided Pauli operators
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

function embed(n::Int,idx,pc::PauliChannel)
    PauliChannel(map(p->(embed(n,idx,p[1]),embed(n,idx,p[2])),pc.paulis), pc.weights)
end

nqubits(pc::PauliChannel) = nqubits(pc.paulis[1][1])

function apply!(state::GeneralizedStabilizer, gate::AbstractPauliChannel; prune_threshold::Union{Nothing, Float64}=nothing)
    if prune_threshold === nothing
        prune_threshold = 1e-14  # Default value
    end
    dict = state.destabweights
    stab = state.stab
    dtype = _dictvaltype(dict)
    tzero = zero(dtype)
    tone = one(dtype)
    newdict = typeof(dict)(tzero)
    for ((dᵢ,dⱼ), χ) in dict # the state
        for ((Pₖ,Pₗ), w) in zip(gate.paulis, gate.weights) # the channel
            phaseₖ, dₖ, dₖˢᵗᵃᵇ = rowdecompose(Pₖ,stab)
            phaseₗ, dₗ, dₗˢᵗᵃᵇ = rowdecompose(Pₗ,stab)
            cₖₗ = (dot(dₖˢᵗᵃᵇ,dᵢ) + dot(dₗˢᵗᵃᵇ,dⱼ))*2
            dᵢ′ = dₖ .⊻ dᵢ
            dⱼ′ = dₗ .⊻ dⱼ
            χ′ = χ * w * (-tone)^cₖₗ * (im)^(-phaseₖ+phaseₗ+4)
            if abs(χ′) >= prune_threshold
                newdict[(dᵢ′,dⱼ′)] = get!(newdict,(dᵢ′,dⱼ′),0)+χ′
            end
        end
    end
    filter!(x -> abs(x[2]) >= prune_threshold, newdict)
    state.destabweights = newdict
    state
end

"""Decompose a Pauli ``P`` in terms of stabilizer and destabilizer rows from a given tableaux.

For given tableaux of rows destabilizer rows ``\\{d_i\\}`` and stabilizer rows ``\\{s_i\\}``,
there are boolean vectors ``b`` and ``c`` such that
``P = i^p \\prod_i d_i^{b_i} \\prod_i s_i^{c_i}``.

By examining the commutation of `P` with the stabilizer and destabilizer generators,
this decomposition can be determined in `Ο(n²)` time.

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
    return p+pauli.phase[], b, c
end

"""A Pauli channel datastructure, mainly for use with [`GeneralizedStabilizer`](@ref).

More convenient to use than [`PauliChannel`](@ref) when you know your Pauli channel is unitary.

```jldoctest
julia> Tgate = UnitaryPauliChannel(
           (I, Z),
           ((1+exp(im*π/4))/2, (1-exp(im*π/4))/2)
       )
A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:
with ϕᵢ | Pᵢ
 0.853553+0.353553im | + _
 0.146447-0.353553im | + Z

julia> PauliChannel(Tgate)
Pauli channel ρ ↦ ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† with the following branches:
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.853553+0.0im | + _ | + _
 0.0+0.353553im | + _ | + Z
 0.0-0.353553im | + Z | + _
 0.146447+0.0im | + Z | + Z
```
"""
struct UnitaryPauliChannel{T,S,P} <: AbstractPauliChannel
    paulis::T
    weights::S
    paulichannel::P # caching the conversion the more general non-unitary type of PauliChannel
    function UnitaryPauliChannel(paulis, weights)
        n = nqubits(paulis[1])
        ws = [w₁*w₂' for w₁ in weights for w₂ in weights]
        ps = [(p₁,p₂) for p₁ in paulis for p₂ in paulis]
        pc = PauliChannel(ps,ws)
        new{typeof(paulis),typeof(weights),typeof(pc)}(paulis,weights,pc)
    end
end

PauliChannel(p::UnitaryPauliChannel) = p.paulichannel

function Base.show(io::IO, pc::UnitaryPauliChannel)
    println(io, "A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:")
    print(io, "with ϕᵢ | Pᵢ")
    for (p, χ) in zip(pc.paulis, pc.weights)
        println(io)
        print(io, " ")
        print(IOContext(io, :compact => true), χ)
        print(io, " | ", p)
    end
end

function embed(n::Int,idx,pc::UnitaryPauliChannel)
    UnitaryPauliChannel(map(p->embed(n,idx,p),pc.paulis), pc.weights)
end

nqubits(pc::UnitaryPauliChannel) = nqubits(pc.paulis[1])

apply!(state::GeneralizedStabilizer, gate::UnitaryPauliChannel; prune_threshold::Union{Nothing, Float64}=nothing) = prune_threshold === nothing ? apply!(state, gate.paulichannel) : apply!(state, gate.paulichannel, prune_threshold)

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> pcT ⊗ P"X"
A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:
with ϕᵢ | Pᵢ
 0.853553+0.353553im | + _X
 0.146447-0.353553im | + ZX
```

"""
function (⊗)(gate::AbstractPauliChannel, Op::PauliOperator)
    new_unitary_channel = typeof(gate)
    ps = typeof(gate.paulichannel)
    new_paulis = pcT.paulis |> collect |> x -> map(y -> y ⊗ Op, x) |> Tuple
    weights = gate.weights
    return UnitaryPauliChannel(new_paulis, weights)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> tHadamard * pcT
A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:
with ϕᵢ | Pᵢ
 0.853553+0.353553im | X₁ ⟼ + Z
Z₁ ⟼ + X
 0.146447-0.353553im | X₁ ⟼ + Z
Z₁ ⟼ - X
```

"""
function (*)(Op::AbstractCliffordOperator, gate::AbstractPauliChannel)
    new_unitary_channel = typeof(gate)
    ps = typeof(gate.paulichannel)
    new_paulis = pcT.paulis |> collect |> x -> map(y -> y * Op, x) |> Tuple
    weights = gate.weights
    return UnitaryPauliChannel(new_paulis, weights)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> sm = GeneralizedStabilizer(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> copy(sm)
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> sm == copy(sm)
true
```

"""
Base.copy(sm::GeneralizedStabilizer) = GeneralizedStabilizer(copy(sm.stab),copy(sm.destabweights))
Base.:(==)(sm₁::GeneralizedStabilizer, sm₂::GeneralizedStabilizer) = sm₁.stab==sm₂.stab && sm₁.destabweights==sm₂.destabweights

##
# Predefined Pauli Channels
##

const pcT = UnitaryPauliChannel(
    (I, Z),
    ((1+exp(im*π/4))/2, (1-exp(im*π/4))/2)
)
