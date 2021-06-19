# Useful States and Operators

```@meta
DocTestSetup = quote
    using QuantumClifford
    using StableRNGs
    rng = StableRNG(42)
end
```

There are numerous frequently used stabilizer states already implemented in this
library.

# Pauli Operators

Single qubit `PauliOperator` is implemented in [`single_z`] and [`single_x`].

```jldoctest
julia> single_z(4,2)
+ _Z__

julia> single_x(4,3)
+ __X_
```

All identity operators use `zero`.

```jldoctest
julia> zero(PauliOperator, 3)
+ ___

julia> zero(P"XYZXYZ")
+ ______
```

Random Pauli operators are implemented as well (with or without a random phase).

```jldoctest rand
julia> using StableRNGs; rng = StableRNG(42);

julia> random_pauli(rng, 4)
+i_ZZ_

julia> random_pauli(rng, 4; nophase=true)
+ ZXZY
```

# Stabilizer States

An all-identity stabilizer can be created with `zero`.

```jldoctest
julia> zero(Stabilizer, 3)
+ ___
+ ___
+ ___

julia> zero(Stabilizer, 2, 3)
+ ___
+ ___

julia> zero(S"XIZ
              YZX")
+ ___
+ ___
```

Diagonal stabilizers in different bases are available as well, through `one`.

```jldoctest
julia> one(Stabilizer, 3)
+ Z__
+ _Z_
+ __Z

julia> one(Stabilizer, 3; basis=:Y)
+ Y__
+ _Y_
+ __Y

julia> one(S"XX
             ZZ")
+ Z_
+ _Z
```

A random stabilizer (or destabilizers or Clifford operators) can be created as well. We use the algorithm described in [bravyi2020hadamard](@cite).

```jldoctest rand
julia> random_stabilizer(rng, 2,5)
+ YZXZZ
- XZYYY
```
## Mixed States

Similarly, one can create a diagonal mixed state.

```jldoctest
julia> one(MixedDestabilizer, 2, 3)
Rank 2 stabilizer
+ X__
+ _X_
━━━━━
+ __X
━━━━━
+ Z__
+ _Z_
━━━━━
+ __Z
```
