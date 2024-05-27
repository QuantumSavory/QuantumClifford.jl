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

# Options for Constructing with MixedDestabilizer

- Given a `Destabilizer` object (which  presumesfull rank), convert it
into a `MixedDestabilizer` object. This allows for the possibility of 
rank deficiency.

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

- Similar to the first option, but with the added capability to
specify the "rank." This rank determines the number of rows
associated with the `Stabilizer` and the number corresponding 
to the logical operators.

```jldoctest mix
julia> MixedDestabilizer(T"ZI IX XX ZZ", 2)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _X
𝒮𝓉𝒶𝒷
+ XX
+ ZZ
```

If the macro string `@T_str` is not convenient, use the normal strings `_T_str`.

```@example mix
julia> MixedDestabilizer(_T_str("ZI IX XX ZZ"), 2)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _X
𝒮𝓉𝒶𝒷
+ XX
+ ZZ
```

- Given a `Stabilizer` object (whichmay have fewer rows than columns)
, perform the necessary canonicalization to determine the 
corresponding destabilizer and logical operators, resulting in a 
complete MixedDestabilizer object.

```jldoctest mix
julia> s = S"-XXX
             +ZZX
             +III";

julia> MixedDestabilizer(s, undoperm=false, reportperm=false)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ _Z_
𝒳ₗ━━━
+ ZZX
𝒮𝓉𝒶𝒷━
+ Y_Y
+ ZXZ
𝒵ₗ━━━
+ Z_Z
```

When `undoperm` is set to `false`, the column permutations are not reversed.
As a result, the qubits may be reindexed according to the permutations
made during the canonicalization process.

```jldoctest mix
julia> MixedDestabilizer(s, undoperm=true, reportperm=false)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ __Z
𝒳ₗ━━━
+ ZXZ
𝒮𝓉𝒶𝒷━
+ YY_
+ ZZX
𝒵ₗ━━━
+ ZZ_
```

When `undoperm` is set to `true`, the column permutations performed during 
canonicalizationare automatically reversed before finalizing the 
`MixedDestabilizer`.

```jldoctest mix
julia> MixedDestabilizer(canonicalize!(s))
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ __Z
𝒳ₗ━━━
+ ZXZ
𝒮𝓉𝒶𝒷━
+ YY_
+ ZZX
𝒵ₗ━━━
+ ZZ_
```

- A low-level constructor that accepts a manually created `Tableau` object. 
Note that the `Tableau` object is not currently public. It serves as the 
underlying data structure for all related objects but does not assume
commutativity or other properties.

```@example mix
julia> MixedDestabilizer(Tableau(Bool[0 0; 0 1; 1 1; 1 0], Bool[1 0; 0 0; 0 0; 1 1]), 2)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z_
+ _X
𝒮𝓉𝒶𝒷
+ XX
+ YZ
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
