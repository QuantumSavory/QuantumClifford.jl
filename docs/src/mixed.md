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
