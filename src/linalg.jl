function permutesystems(c::CliffordOperator,p) # TODO this is a slow stupid implementation
    CliffordOperator(Tableau([tab(c)[i][p] for i in 1:2*nqubits(c)][vcat(p,p.+nqubits(c))]))
end

@deprecate permute(c::CliffordOperator,p) permutesystems(c,p)

function permutesystems!(s::Tableau, perm::AbstractVector)
    for r in 1:size(s,1)
        s[r] = s[r][perm] # TODO make a local temporary buffer row instead of constantly allocating new rows
    end
    s
end

function permutesystems!(s::AbstractStabilizer, perm::AbstractVector)
    permutesystems!(tab(s), perm)
    s
end

import Base: permute!
@deprecate permute!(s::Tableau, perm::AbstractVector) permutesystems!(s, perm)
@deprecate permute!(s::AbstractStabilizer, perm::AbstractVector) permutesystems!(s, perm)

# TODO upstream to QuantumInterface for (state::Any, perm)
permutesystems(s::AbstractStabilizer, perm) = permutesystems!(copy(s), perm)

"""
$TYPEDSIGNATURES

Inverse of a `CliffordOperator`

```jldoctest
julia> inv(CliffordOperator(sCNOT))
X₁ ⟼ + XX
X₂ ⟼ + _X
Z₁ ⟼ + Z_
Z₂ ⟼ + ZZ

julia> inv(CliffordOperator(sCNOT(2, 1), 2))
X₁ ⟼ + X_
X₂ ⟼ + XX
Z₁ ⟼ + ZZ
Z₂ ⟼ + _Z

julia> inv(CliffordOperator(tHadamard))
X₁ ⟼ + Z
Z₁ ⟼ + X
```
"""
function LinearAlgebra.inv(c::CliffordOperator; phases=true)
    ci = zero(c)
    n = nqubits(c)
    # TODO this transpose can be much faster with proper SIMDing
    for i in 1:n
        for j in 1:n
            ci.tab[i,j] = c.tab[n+j,i][2], c.tab[j,i][2]
            ci.tab[n+i,j] = c.tab[n+j,i][1], c.tab[j,i][1]
        end
    end
    if phases
        ci*c*ci # TODO perform this inplace as in Stim https://github.com/quantumlib/Stim/blob/e51ea66d213b25920e72c08e53266ec56fd14db4/src/stim/stabilizers/tableau.cc#L383
    else
        ci
    end
end

"""The inner product of two Stabilizers.

Based on [garcia2012efficient](@cite).

```jldoctest
julia> using LinearAlgebra

julia> dot(S"Z", S"Z")
1.0

julia> dot(S"Z", S"Y")
0.7071067811865476
```

See also: [`logdot`](@ref)"""
function LinearAlgebra.dot(s1::AbstractStabilizer, s2::AbstractStabilizer)
    ld = logdot(s1,s2)
    if isnothing(ld)
        return 0.0
    else
        return 2.0^(-ld/2)
    end
end

"""Logarithm of the inner product between to Stabilizer states.

If the result is `nothing`, the dot inner product is zero.
Otherwise the inner product is `2^(-logdot/2)`.

The actual inner product can be computed with `LinearAlgebra.dot`.

Based on [garcia2012efficient](@cite)."""
function logdot(s1::AbstractStabilizer, s2::AbstractStabilizer) # TODO verify rank # TODO this is currently very inefficient as we discard the destabilizers and then recreate them
    logdot(stabilizerview(s1),stabilizerview(s2))
end

function logdot(s1::Stabilizer, s2::Stabilizer)
    if nqubits(s1)!=length(s1) || nqubits(s2)!=length(s2) # TODO implement this
        throw(DomainError("Only pure (not mixed) states are supported when calculating inner product."))
    end
    if nqubits(s1)!=nqubits(s2)
        throw(DimensionMismatch("Inner product can be calculated only between states with the same number of qubits."))
    end
    c1_inv = inv(CliffordOperator(tab(MixedDestabilizer(copy(s1)))))
    s2_prime = canonicalize!(c1_inv*s2)
    canonicalize!(s2_prime)
    k = 0
    for row in eachindex(s2_prime)
        if any(i->s2_prime[row,i][1], 1:nqubits(s2_prime)) # X or Y
            k += 1
        else
            if !iszero(tab(s2_prime).phases[row])
                return nothing
            end
        end
    end
    return k
end

"""
    fidelity(state1::MixedDestabilizer, state2::MixedDestabilizer)

Compute the fidelity between stabilizer states, following the
`QuantumInterface.fidelity` convention,

```math
F(\\rho, \\sigma) =
\\operatorname{Tr}\\sqrt{\\sqrt{\\rho}\\,\\sigma\\sqrt{\\rho}}.
```

For two pure states this returns their overlap by delegating to
`LinearAlgebra.dot`. If exactly one state is pure, the result is the square root
of that pure state's projector expectation in the mixed state. The order of
the arguments does not matter. Inputs have the same rank and representation
requirements as [`expect`](@ref).

General mixed/mixed fidelity is not currently implemented and raises a
`DomainError`. In particular, it is not equal in general to the square root of
the Hilbert--Schmidt expectation returned by [`expect`](@ref).

# Examples

```jldoctest
julia> fidelity(MixedDestabilizer(S"Z"), maximally_mixed(1))
0.7071067811865476

julia> fidelity(MixedDestabilizer(S"Z"), MixedDestabilizer(S"-Z"))
0.0
```

See also: [`expect`](@ref), [`logdot`](@ref).
"""
function fidelity(
    state1::MixedDestabilizer,
    state2::MixedDestabilizer,
)
    state1_qubits = nqubits(state1)
    state2_qubits = nqubits(state2)
    state1_qubits == state2_qubits || throw(DimensionMismatch(
        lazy"Expected equal qubit counts; got $state1_qubits and $state2_qubits.",
    ))
    state1_rank = rank(state1)
    state2_rank = rank(state2)
    state1_is_pure = state1_rank == state1_qubits
    state2_is_pure = state2_rank == state2_qubits

    if state1_is_pure && state2_is_pure
        return dot(state1, state2)
    elseif !state1_is_pure && !state2_is_pure
        throw(DomainError(
            (state1_rank, state2_rank),
            lazy"Fidelity between two rank-deficient stabilizer states is not supported.",
        ))
    end

    pure_state, mixed_state = state1_is_pure ?
        (state1, state2) : (state2, state1)
    exponent = _stabilizer_expect_log2(nothing, pure_state, mixed_state)
    isnothing(exponent) ? 0.0 : exp2(exponent / 2)
end

LinearAlgebra.rank(s::Stabilizer)   = throw(BadDataStructure("Using a `Stabilizer` type does not permit automatic tracking of the rank. Use `length`, `trusted_rank`, the `MixedDestabilizer` type, or track the rank manually.",
                                            :rank, :Stabilizer))
LinearAlgebra.rank(s::Destabilizer) = throw(BadDataStructure("Using a `Destabilizer` type does not permit automatic tracking of the rank. Use `length`, `trusted_rank`, the `MixedDestabilizer` type, or track the rank manually.",
                                            :rank, :Destabilizer))
LinearAlgebra.rank(s::MixedStabilizer) = s.rank
LinearAlgebra.rank(s::MixedDestabilizer) = s.rank

"""A "trusted" `rank` which returns `rank(state)` for `Mixed[De]Stabilizer` and `length(state)` for `[De]Stabilizer`."""
function trusted_rank end
trusted_rank(s::Stabilizer) = length(s)
trusted_rank(s::Destabilizer) = length(s)
trusted_rank(s::MixedStabilizer) = LinearAlgebra.rank(s)
trusted_rank(s::MixedDestabilizer) = LinearAlgebra.rank(s)

"""Tensor product between operators or tableaux.

Tensor product between CiffordOperators:

```jldoctest
julia> tensor(CliffordOperator(sCNOT), CliffordOperator(sCNOT))
X₁ ⟼ + XX__
X₂ ⟼ + _X__
X₃ ⟼ + __XX
X₄ ⟼ + ___X
Z₁ ⟼ + Z___
Z₂ ⟼ + ZZ__
Z₃ ⟼ + __Z_
Z₄ ⟼ + __ZZ
```

Tensor product between PauliOperators:

```jldoctest
julia> tensor(P"-IXYZ", P"iIXYZ")
-i_XYZ_XYZ
```

Tensor product between Tableaux:

```jldoctest
julia> s = S"-XX
             +ZZ";

julia> tensor(s, s, s)
- XX____
+ ZZ____
- __XX__
+ __ZZ__
- ____XX
+ ____ZZ

julia> s = S"+XZI
             -IZI";

julia> tensor(s, s)
+ XZ____
- _Z____
+ ___XZ_
- ____Z_
```

See also [`tensor_pow`](@ref)."""
function tensor end

tensor(p::PauliOperator, ps::PauliOperator...) = foldl(⊗, ps; init=p)

function tensor(ops::AbstractStabilizer...) # TODO optimize by pre-allocating one large tableau instead of the current quadratic fold
    ct = promote_type(map(typeof, ops)...)
    conv_ops = map(x -> convert(ct, x), ops)
    return foldl(⊗, conv_ops)
end

"""Repeated tensor product of an operators or a tableau.

For `CliffordOperator`:

```jldoctest
julia> tensor_pow(CliffordOperator(sHadamard), 3)
X₁ ⟼ + Z__
X₂ ⟼ + _Z_
X₃ ⟼ + __Z
Z₁ ⟼ + X__
Z₂ ⟼ + _X_
Z₃ ⟼ + __X
```

For `PauliOperator`:

```jldoctest
julia> tensor_pow(P"IXYZ", 2)
+ _XYZ_XYZ
```

For `Tableaux`:

```jldoctest
julia> tensor_pow(S"Z", 4)
+ Z___
+ _Z__
+ __Z_
+ ___Z

julia> s = S"+XZI
             +IZI";

julia> tensor_pow(s, 3)
+ XZ_______
+ _Z_______
+ ___XZ____
+ ____Z____
+ ______XZ_
+ _______Z_
```

See also [`tensor`](@ref).
"""
function tensor_pow(op::Union{<:AbstractStabilizer,<:AbstractCliffordOperator},power)
    if power==1
        return op
    else
        return tensor((op for i in 1:power)...)
    end
end

function tensor(ops::Stabilizer...)
    length(ops)==1 && return ops[1]
    ntot = sum(nqubits, ops) # TODO why is this allocating (at least in 1.11)
    rtot = sum(length, ops)  # TODO why is this allocating (at least in 1.11)
    tab = zero(Stabilizer, rtot, ntot)
    last_row = 0
    last_col = 0
    for op in ops
        _, last_row, last_col = puttableau!(tab, op, last_row, last_col)
    end
    tab
end

function tensor(ops::MixedDestabilizer...)
    length(ops)==1 && return ops[1]
    ntot = sum(nqubits, ops)
    rtot = sum(LinearAlgebra.rank, ops)
    tab = zero(Tableau, 2*ntot, ntot)
    last_svrow = ntot
    last_dvrow = 0
    last_lxrow = rtot
    last_lzrow = ntot+rtot
    last_col = 0
    for op in ops
        _, last_svrow, _        = puttableau!(tab,   stabilizerview(op), last_svrow, last_col)
        _, last_dvrow, _        = puttableau!(tab, destabilizerview(op), last_dvrow, last_col)
        _, last_lxrow, _        = puttableau!(tab,     logicalxview(op), last_lxrow, last_col)
        _, last_lzrow, last_col = puttableau!(tab,     logicalzview(op), last_lzrow, last_col)
    end
    MixedDestabilizer(tab, rtot)
end

tensor(ops::AbstractCliffordOperator...) = ⊗(CliffordOperator.(ops)...)

function tensor(ops::CliffordOperator...) # TODO implement \otimes for Destabilizer and use it here
    length(ops)==1 && return ops[1]
    ntot = sum(nqubits, ops)
    tab = zero(Tableau, 2*ntot, ntot)
    last_zrow = ntot
    last_xrow = 0
    for op in ops
        t = QuantumClifford.tab(op)
        _, last_zrow, _ = puttableau!(tab, (@view t[end÷2+1:end]), last_zrow, last_xrow)
        _, last_xrow, _ = puttableau!(tab, (@view t[1:end÷2]), last_xrow, last_xrow)
    end
    CliffordOperator(tab)
end
