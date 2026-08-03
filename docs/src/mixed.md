# [Mixed Stabilizer States](@id Mixed-Stabilizer-States)

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

The [`Stabilizer`](@ref) and [`Destabilizer`](@ref) have some support for mixed
states (by being initialized with an incomplete list of stabilizer generators),
but for most purposes one would use the `Mixed*` data structures.

Mixed stabilizer states are implemented with [`MixedStabilizer`](@ref) and
[`MixedDestabilizer`](@ref), the latter of which is the preferred data structure
for most tasks as it is much faster by virtue of tracking the destabilizer
generators.

```jldoctest mix
julia> s = S"XXX
             IZZ";

julia> Destabilizer(s)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ _X_
𝒮𝓉𝒶𝒷━
+ XXX
+ _ZZ
```

Unlike `Destabilizer`, `MixedDestabilizer` also tracks the logical
operation generators.

```jldoctest mix
julia> m = MixedDestabilizer(s)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ _X_
𝒳ₗ━━━
+ _XX
𝒮𝓉𝒶𝒷━
+ XXX
+ _ZZ
𝒵ₗ━━━
+ Z_Z

julia> stabilizerview(m)
+ XXX
+ _ZZ

julia> destabilizerview(m)
+ Z__
+ _X_

julia> logicalxview(m)
+ _XX

julia> logicalzview(m)
+ Z_Z
```

# Density Operators, Expectations, and Fidelity

An `n`-qubit stabilizer tableau with `r` independent generators represents the
normalized density operator

```math
\rho_S = 2^{-n}\sum_{g\in S}g = \frac{P_S}{2^{n-r}},
```

where ``S`` is the signed stabilizer group and ``P_S`` projects onto its common
``+1`` eigenspace. A rank-`n` tableau is pure. A rank-deficient tableau is the
maximally mixed state on a stabilized subspace of dimension ``2^{n-r}``.

[`expect`](@ref) accepts a stabilizer state as the operator as well as a
[`PauliOperator`](@ref). For two state arguments it returns their
Hilbert--Schmidt expectation ``\operatorname{Tr}(\rho_A\rho_B)``. When the
operator state is pure, this is the probability of its projector.

```jldoctest mix
julia> mixed = maximally_mixed(1);

julia> expect(S"Z", mixed)
0.5

julia> expect(mixed, mixed) # purity of I/2
0.5
```

The indexed form embeds the operator at an ordered list of subsystems, leaving
the identity on all other qubits. The order is significant: operator qubit 1
maps to `indices[1]`, operator qubit 2 maps to `indices[2]`, and so on. This
computes local projector probabilities without first permuting or tracing out
the state.

```jldoctest mix
julia> expect([1, 3], bell(), ghz(3))
0.5
```

[`fidelity`](@ref) uses the root-fidelity convention. For two pure states it
returns their overlap and delegates to `LinearAlgebra.dot` when their tableau
presentations satisfy that function's input contract. If exactly one state is
pure, its square is the corresponding pure-state projector expectation.

```jldoctest mix
julia> fidelity(S"Z", mixed)
0.7071067811865476
```

General mixed/mixed Uhlmann fidelity needs more information than the
Hilbert--Schmidt expectation and is not yet implemented. Calling `fidelity`
with two rank-deficient tableaux raises a `DomainError` rather than returning a
different quantity under the fidelity name.

# Gottesman Canonicalization

To obtain the logical operators we perform a different type of canonicalization,
described in Gottesman's thesis and implemented in [`canonicalize_gott!`](@ref).
Unlike [`canonicalize!`](@ref) which uses only row operations,
`canonicalize_gott!` performs column swaps as well. `MixedDestabilizer` undoes
those swaps by default when instantiated, but that behavior can be turned off,
if you prefer to work with the canonicalized tableau.

```jldoctest mix
julia> s = S"XXX
             ZIZ";

julia> MixedDestabilizer(s)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ __X
𝒳ₗ━━━
+ _X_
𝒮𝓉𝒶𝒷━
+ XXX
+ Z_Z
𝒵ₗ━━━
+ ZZ_

julia> MixedDestabilizer(s; undoperm=false)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ _X_
𝒳ₗ━━━
+ __X
𝒮𝓉𝒶𝒷━
+ XXX
+ ZZ_
𝒵ₗ━━━
+ Z_Z
```

`Destabilizer` and `MixedStabilizer` do not use any column swaps on
instantiation as they do not track the logical operators.
