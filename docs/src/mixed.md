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
+ Z__
+ _X_
━━━━━
+ XXX
+ _ZZ
```

Unlike `Destabilizer`, `MixedDestabilizer` also tracks the logical
operation generators.

```jldoctest mix
julia> m = MixedDestabilizer(s)
Rank 2 stabilizer
+ Z__
+ _X_
━━━━━
+ _XX
━━━━━
+ XXX
+ _ZZ
━━━━━
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
Rank 2 stabilizer
+ Z__
+ __X
━━━━━
+ _X_
━━━━━
+ XXX
+ Z_Z
━━━━━
+ ZZ_

julia> MixedDestabilizer(s; undoperm=false)
Rank 2 stabilizer
+ Z__
+ _X_
━━━━━
+ __X
━━━━━
+ XXX
+ ZZ_
━━━━━
+ Z_Z
```

`Destabilizer` and `MixedStabilizer` do not use any column swaps on
instantiation as they do not track the logical operators.
