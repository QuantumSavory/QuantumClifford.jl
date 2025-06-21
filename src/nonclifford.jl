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

julia> pcT
A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:
with ϕᵢ | Pᵢ
 0.853553+0.353553im | + _
 0.146447-0.353553im | + Z

julia> apply!(GeneralizedStabilizer(S"-X"), pcT)
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
```

See also: [`PauliChannel`](@ref)
"""
mutable struct GeneralizedStabilizer{T,F} <: AbstractQCState
    stab::T
    destabweights::DefaultDict{Tuple{BitVector, BitVector}, F, F}
end

function GeneralizedStabilizer(state)
    n = nqubits(state)
    md = MixedDestabilizer(state)
    rank(md)==n || throw(ArgumentError(lazy"""Attempting to convert a `Stabilizer`-like object to `GeneralizedStabilizer` object failed, because the initial state does not represent a pure state. Currently only pure states can be used to initialize a `GeneralizedStabilizer` mixture of stabilizer states."""))
    GeneralizedStabilizer(md, DefaultDict(0.0im, (falses(n),falses(n))=>1.0+0.0im)) # TODO maybe it should default to Destabilizer, not MixedDestabilizer
end

GeneralizedStabilizer(s::GeneralizedStabilizer) = s
Base.copy(sm::GeneralizedStabilizer) = GeneralizedStabilizer(copy(sm.stab),copy(sm.destabweights))
Base.:(==)(sm₁::GeneralizedStabilizer, sm₂::GeneralizedStabilizer) = sm₁.stab==sm₂.stab && sm₁.destabweights==sm₂.destabweights

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

"""
Apply a Clifford gate to a generalized stabilizer state, i.e. a weighted sum of stabilizer states.

```jldoctest
julia> sm = GeneralizedStabilizer(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> apply!(sm, CliffordOperator(tHadamard))
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ X
𝒮𝓉𝒶𝒷
- Z
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _
```

See also: [`GeneralizedStabilizer`](@ref)
"""
function apply!(state::GeneralizedStabilizer, gate::AbstractCliffordOperator)
    apply!(state.stab, gate)
    state
end

"""$(TYPEDSIGNATURES)

Expectation value for the [PauliOperator](@ref) observable given the [`GeneralizedStabilizer`](@ref) state `s`.

```jldoctest genstab
julia> sm = GeneralizedStabilizer(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> apply!(sm, pcT)
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

julia> χ′ = expect(P"-X", sm)
0.7071067811865475 + 0.0im

julia> prob = (real(χ′)+1)/2
0.8535533905932737
```

"""
function expect(p::PauliOperator, s::GeneralizedStabilizer) # TODO optimize
    χ′ = zero(valtype(s.destabweights))
    phase, b, c = rowdecompose(p, s.stab)
    for ((dᵢ,dⱼ), χ) in s.destabweights
        _allthreesumtozero(dᵢ,dⱼ,b) || continue
        χ′ += χ * (-1)^(dᵢ'*c)
    end
    return (im)^(phase) * χ′
end

"""Same as `all(==(0), (a.+b.+c) .% 2)`"""
function _allthreesumtozero(a,b,c)
    n = length(a)
    @inbounds @simd for i in 1:n
        odd = (a[i]+b[i]+c[i]) & 1
        if odd != 0
            return false
        end
    end
    true
end

"""Compute the trace of a [`GeneralizedStabilizer`](@ref) state.

```jldoctest trace
julia> using QuantumClifford; using LinearAlgebra;

julia> sm = GeneralizedStabilizer(S"-X");

julia> apply!(sm, pcT)
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

julia> tr(sm)
1.0 + 0.0im
```
"""
function LinearAlgebra.tr(sm::GeneralizedStabilizer)
    trace_χ′ = sum(χ for ((P_i, P_j), χ) in sm.destabweights if P_i == P_j; init=0)
    return trace_χ′
end

"""Returns the updated `GeneralizedStabilizer` state sm′ = (χ′, B(S′, D′)),
where (S′, D′) is derived from (S, D) through the traditional stabilizer update,
and χ′ is the updated density matrix after measurement. Note: Λ(χ′) ≤ Λ(χ).
"""
function _projectrand_notnorm(sm::GeneralizedStabilizer, p::PauliOperator, res::Int)
    dict = sm.destabweights
    dtype = valtype(dict)
    tzero = zero(dtype)
    tone = one(dtype)
    stab = sm.stab
    newdict = typeof(dict)(tzero)
    phase, b, c = rowdecompose(p, stab)
    new_stab = copy(stab)
    s_view = stabilizerview(new_stab)
    d_view = destabilizerview(new_stab)
    n = nqubits(new_stab)
    id_op = zero(PauliOperator, n)
    sign = res == 0 ? 1 : -1

    # Implementation of the in-place Pauli measurement quantum operation (Algorithm 2)
    # on a generalized stabilizer by Ted Yoder (Page 8) from [Yoder2012AGO](@cite).
    if all(iszero, b)
        # (Eq. 14-17)
        for ((dᵢ, dⱼ), χ) in dict
            cond_i = im^phase * (-tone)^(dot(dᵢ, c)) == sign # (Eq. 16)
            cond_j = im^phase * (-tone)^(dot(dⱼ, c)) == sign # (Eq. 16)
            if cond_i && cond_j
                newdict[(dᵢ, dⱼ)] += χ
            end
        end
        sm.destabweights = newdict
        return sm, res
    else
        # (Eq. 18-26)
        k_pos = findfirst(!iszero, b)
        # get the k-th stabilizer generator
        sk = s_view[k_pos]
        # update stabilizer generators
        for j in eachindex(b)
            if b[j]
                s_view[j] *= sk
            end
        end
        # update destabilizer generators
        for j in eachindex(c)
            if c[j] && j != k_pos  # cj = 1 and j ≠ k_pos
                d_view[j] *= sk
            end
        end
        d_view[k_pos] = id_op
        d_view[k_pos] = sk
        # replace sk with sign*M
        s_view[k_pos] = sign * p
        # update the χ matrix
        k = falses(n)
        k[k_pos] = true
        for ((dᵢ, dⱼ), χ) in dict
            x, y = dᵢ, dⱼ
            q = one(dtype)
            if dot(dᵢ, k) == 1
                q *= im^phase * (-tone)^dot(dᵢ, c) * sign
                x = dᵢ .⊻ b
            end
            if dot(dⱼ, k) == 1
                q *= conj(im^phase) * (-tone)^dot(dⱼ, c) * sign
                y = dⱼ .⊻ b
            end
            newdict[(x, y)] += χ * q / 2
        end
        sm.destabweights = newdict
        sm.stab = new_stab
        return sm, res
    end
end

"""$(TYPEDSIGNATURES)

Performs a randomized projection of the state represented by the [`GeneralizedStabilizer`](@ref) `sm`,
based on the measurement of a [PauliOperator](@ref) `p`.

Unlike in the case of stabilizer states, the expectation value χ′ of a Pauli operator
with respect to these more general states can be any real number between -1 and 1.
The expectation value can be calculated with `expect(p, sm)`.

```math
\\chi' = \\langle p \\rangle = \\text{expect}(p, sm)
```

To convert χ′ into a probability of projecting on the +1 eigenvalue branch:

```math
\\text{probability}_{1} = \\frac{\\text{real}(\\chi') + 1}{2}
```

!!! note
    Because the possible measurement results are themselves not stabilizer states anymore,
    we can not use the `project!` API, which assumes a stabilizer tableau and reports detailed
    information about whether the tableau and measurement commute or anticommute.

```jldoctest genstab
julia> sm = GeneralizedStabilizer(S"-X");

julia> apply!(sm, pcT)
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

julia> χ′ = expect(P"-X", sm)
0.7071067811865475 + 0.0im

julia> prob₁ = (real(χ′)+1)/2
0.8535533905932737

julia> QuantumClifford._projectrand_notnorm(copy(sm), P"X", 0)[1]
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.146447+0.0im | + Z | + Z

julia> QuantumClifford._projectrand_notnorm(copy(sm), P"X", 1)[1]
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.853553+0.0im | + _ | + _
```

See also: [`expect`](@ref)
"""
function projectrand!(sm::GeneralizedStabilizer, p::PauliOperator)
    # Compute expectation value
    exp_val = expect(p, sm)
    prob_plus = (real(exp_val) + 1) / 2
    # Randomly choose outcome
    res = rand() < prob_plus ? 0 : 1
    # Apply the corresponding projection
    sm, _ = _projectrand_notnorm(sm, p, res)
    # Normalize the state
    trace = tr(sm)
    for ((dᵢ, dⱼ), χ) in sm.destabweights
        sm.destabweights[(dᵢ, dⱼ)] = χ / trace
    end
    res = (res == 0) ? 0x0 : 0x2
    return sm, res
end

function project!(s::GeneralizedStabilizer, p::PauliOperator)
    throw(MethodError(project!, (s, p)))
end

nqubits(sm::GeneralizedStabilizer) = nqubits(sm.stab)

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
    for (i, (di,dj)) in enumerate(pc.paulis)
        χ = pc.weights[i]
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

"""Applies a (potentially non-unitary) Pauli channel to a generalized stabilizer.

See also: [`GeneralizedStabilizer`](@ref), [`PauliChannel`](@ref), [`UnitaryPauliChannel`](@ref)
"""
function apply!(state::GeneralizedStabilizer, gate::AbstractPauliChannel; prune_threshold=1e-10)
    nqubits(state) == nqubits(gate) || throw(DimensionMismatch(lazy"GeneralizedStabilizer has $(nqubits(state)) qubits, but PauliChannel has $(nqubits(gate)). Use `embed` to create an appropriately padded PauliChannel."))
    dict = state.destabweights
    stab = state.stab
    dtype = valtype(dict)
    tzero = zero(dtype)
    tone = one(dtype)
    newdict = typeof(dict)(tzero)
    for ((dᵢ,dⱼ), χ) in dict # the state
        for (i, (Pₗ,Pᵣ)) in enumerate(gate.paulis) # the channel
            w = gate.weights[i]
            phaseₗ, dₗ, dₗˢᵗᵃᵇ = rowdecompose(Pₗ,stab)
            phaseᵣ, dᵣ, dᵣˢᵗᵃᵇ = rowdecompose(Pᵣ,stab)
            c = (dot(dₗˢᵗᵃᵇ,dᵢ) + dot(dᵣˢᵗᵃᵇ,dⱼ))*2
            dᵢ′ = dₗ .⊻ dᵢ
            dⱼ′ = dᵣ .⊻ dⱼ
            χ′ = χ * w * (-tone)^c * (im)^(-phaseₗ+phaseᵣ+4)
            newdict[(dᵢ′,dⱼ′)] += χ′
        end
    end
    filter!(x -> abs(x[2]) > prune_threshold, newdict)
    state.destabweights = newdict
    state
end

"""
$(TYPEDSIGNATURES)

Tensor product of [`GeneralizedStabilizer`](@ref) states.

# Stabilizer state

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

# Arbitrary state

```jldoctest
julia> using LinearAlgebra; # hide

julia> sm = GeneralizedStabilizer(ghz(2))
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _X
𝒮𝓉𝒶𝒷
+ XX
+ ZZ
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + __ | + __

julia> apply!(sm, embed(2, 2, pcT))
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _X
𝒮𝓉𝒶𝒷
+ XX
+ ZZ
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.853553+0.0im | + __ | + __
 0.0+0.353553im | + __ | + Z_
 0.0-0.353553im | + Z_ | + __
 0.146447+0.0im | + Z_ | + Z_

julia> newsm = sm ⊗ sm
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z___
+ _X__
+ __Z_
+ ___X
𝒮𝓉𝒶𝒷━━
+ XX__
+ ZZ__
+ __XX
+ __ZZ
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.0-0.301777im | + Z___ | + ____
 -0.125+0.0im | + Z_Z_ | + ____
 0.125+0.0im | + Z___ | + Z___
 0.728553+0.0im | + ____ | + ____
 0.0-0.0517767im | + Z_Z_ | + Z___
 0.0-0.301777im | + __Z_ | + ____
 0.0+0.301777im | + ____ | + Z___
 0.125+0.0im | + __Z_ | + Z___
 0.125+0.0im | + Z___ | + __Z_
 0.0-0.0517767im | + Z_Z_ | + __Z_
 0.0+0.0517767im | + Z___ | + Z_Z_
 0.0+0.301777im | + ____ | + __Z_
 0.0214466+0.0im | + Z_Z_ | + Z_Z_
 0.125+0.0im | + __Z_ | + __Z_
 -0.125+0.0im | + ____ | + Z_Z_
 0.0+0.0517767im | + __Z_ | + Z_Z_

julia> real(tr(newsm))
1.0
```
"""
function (⊗)(state₁::GeneralizedStabilizer, state₂::GeneralizedStabilizer)
    dict₁ = state₁.destabweights
    dict₂ = state₂.destabweights
    dtype = valtype(dict₁)
    tzero = zero(dtype)
    newdict = DefaultDict{Tuple{BitVector,BitVector},dtype}(tzero)
    newstab = state₁.stab ⊗ state₂.stab
    for ((d1_i, d1_j), χ) in dict₁ # χ = ϕᵢⱼ for state₁
        for ((d2_i, d2_j), χ′) in dict₂ # χ′ = ϕₖₗ for state₂
            # Combine the Pauli operators via tensor product: Pᵢ ⊗ Pₖ
            # and Pⱼ ⊗ Pₗ. vcat implements P₁ ⊗ P₂ as bitwise concatenation.
            new_key_i = vcat(d1_i, d2_i)
            new_key_j = vcat(d1_j, d2_j)
            # The new coefficient is ϕᵢⱼ * ϕₖₗ because:
            # ∑ϕᵢⱼ*ϕₖₗ(Pᵢρ₁Pⱼ†) ⊗ (Pₖρ₂Pₗ†) = ∑ϕᵢⱼ*ϕₖₗ(Pᵢ ⊗ Pₖ)(ρ₁ ⊗ ρ₂)(Pⱼ ⊗ Pₗ)† and
            # thus the combined weight is the product of the individual weights.
            newdict[(new_key_i, new_key_j)] += χ * χ′
        end
    end
    return GeneralizedStabilizer(newstab, newdict)
end

"""Tensor product between [`GeneralizedStabilizer`](@ref) and [`Stabilizer`](@ref).

```jldoctest
julia> using LinearAlgebra; # hide

julia> sm = GeneralizedStabilizer(ghz(2))
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _X
𝒮𝓉𝒶𝒷
+ XX
+ ZZ
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + __ | + __

julia> apply!(sm, embed(2, 2, pcT))
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _X
𝒮𝓉𝒶𝒷
+ XX
+ ZZ
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.853553+0.0im | + __ | + __
 0.0+0.353553im | + __ | + Z_
 0.0-0.353553im | + Z_ | + __
 0.146447+0.0im | + Z_ | + Z_

julia> s = ghz(2)
+ XX
+ ZZ

julia> newsm = sm ⊗ s
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z___
+ _X__
+ __Z_
+ ___X
𝒮𝓉𝒶𝒷━━
+ XX__
+ ZZ__
+ __XX
+ __ZZ
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.0+0.353553im | + ____ | + Z___
 0.0-0.353553im | + Z___ | + ____
 0.146447+0.0im | + Z___ | + Z___
 0.853553+0.0im | + ____ | + ____

julia> real(tr(newsm))
1.0
```
"""
tensor(ops::Union{AbstractStabilizer,GeneralizedStabilizer}...) = tensor(GeneralizedStabilizer.(ops)...)

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
UnitaryPauliChannel(p::UnitaryPauliChannel) = p
UnitaryPauliChannel(p::PauliOperator) = UnitaryPauliChannel([p], [1.0+0.0im])
PauliChannel(p::PauliOperator) = UnitaryPauliChannel(p).paulichannel
PauliChannel(p::UnitaryPauliChannel) = p.paulichannel
Base.copy(p::UnitaryPauliChannel) = UnitaryPauliChannel(map(copy, p.paulis), map(copy, p.weights))
Base.:(==)(p₁::UnitaryPauliChannel, p₂::UnitaryPauliChannel) = p₁.paulis==p₂.paulis && p₁.weights==p₂.weights

function Base.show(io::IO, pc::UnitaryPauliChannel)
    println(io, "A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:")
    print(io, "with ϕᵢ | Pᵢ")
    for (i, p) in enumerate(pc.paulis)
        χ = pc.weights[i]
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

apply!(state::GeneralizedStabilizer, gate::UnitaryPauliChannel; prune_threshold=1e-10) = apply!(state, gate.paulichannel; prune_threshold)

function tensor(pcs::UnitaryPauliChannel...)
    newpaulis = [tensor(ps...) for ps in Iterators.product([pc.paulis for pc in pcs]...)]
    newweights = [prod(ws) for ws in Iterators.product([pc.weights for pc in pcs]...)]
    return UnitaryPauliChannel(newpaulis, newweights)
end

"""
Tensor product between [`UnitaryPauliChannel`](@ref) and [`PauliOperator`](@ref).

```jldoctest
julia> pcT ⊗ P"X"
A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:
with ϕᵢ | Pᵢ
 0.853553+0.353553im | + _X
 0.146447-0.353553im | + ZX

julia> pcT ⊗ pcT ⊗ P"X"
A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:
with ϕᵢ | Pᵢ
 0.603553+0.603553im | + __X
 0.25-0.25im | + Z_X
 0.25-0.25im | + _ZX
 -0.103553-0.103553im | + ZZX
```
"""
tensor(pcs::Union{UnitaryPauliChannel,PauliOperator}...) = tensor(UnitaryPauliChannel.(pcs)...)

"""
Apply a [`UnitaryPauliChannel`](@ref) to a [`GeneralizedStabilizer`](@ref) state.

```jldoctest
julia> sm = GeneralizedStabilizer(S"-X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
- X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> pcT*sm
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
```
"""
Base.:(*)(pc::UnitaryPauliChannel, sm::GeneralizedStabilizer) = apply!(sm, pc.paulichannel)

"""
Calculates the number of non-zero elements in the density matrix `χ`
of a [`GeneralizedStabilizer`](@ref), representing the inverse sparsity
of `χ`. It provides a measure of the state's complexity, with bounds
`Λ(χ) ≤ 4ⁿ`.

```jldoctest heuristic
julia> using QuantumClifford: invsparsity; # hide

julia> sm = GeneralizedStabilizer(S"X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
+ X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> apply!(sm, pcT) |> invsparsity
4
```

Similarly, it calculates the number of non-zero elements in the density
matrix `ϕᵢⱼ`​ of a PauliChannel, providing a measure of the channel
complexity.

```jldoctest heuristic
julia> invsparsity(pcT)
4
```

See also: [`GeneralizedStabilizer`](@ref)
"""
function invsparsity end

invsparsity(sm::GeneralizedStabilizer) = count(!iszero, values(sm.destabweights::DefaultDict{Tuple{BitVector, BitVector}, ComplexF64, ComplexF64}))
invsparsity(gate::AbstractPauliChannel) = count(!iszero, values(gate.paulichannel.weights::Vector{ComplexF64}))

##
# Predefined Pauli Channels
##

const pcT = UnitaryPauliChannel(
    (I, Z),
    ((1+exp(im*π/4))/2, (1-exp(im*π/4))/2)
)
