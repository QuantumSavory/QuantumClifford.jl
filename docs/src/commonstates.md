# [Useful States and Operators](@id Useful-States-and-Operators)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using StableRNGs
    rng = StableRNG(42)
end
```
## States

Stabilizer states can be represented with the [`Stabilizer`](@ref), [`Destabilizer`](@ref), [`MixedStabilizer`](@ref), and [`MixedDestabilizer`](@ref) tableau data structures. You probably want to use [`MixedDestabilizer`](@ref) which supports the widest set of operations.

Moreover, a [`MixedDestabilizer`](@ref) can be stored inside a [`Register`](@ref) together with a set of classical bits in which measurement results can be written.

Below are convenience constructors for common types of states and operators,
already implemented in this library.

## Pauli Operators

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

## Stabilizer States

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

## Enumerating all Clifford Operations

The algorithm from [koenig2014efficiently](@cite) can be used to enumerate all Clifford operations on a given number of qubits through [`enumerate_cliffords`](@ref).
Or one can use [random_clifford](@ref), [random_stabilizer](@ref) to directly sample from that set.

```jldoctest
julia> length(enumerate_cliffords(1))
6

julia> length(enumerate_cliffords(2))
720
```

To also enumerate possible phases, you can use [`enumerate_phases`](@ref).

```jldoctest
julia> length(collect(enumerate_phases(tCNOT)))
16

julia> length(collect(enumerate_phases(enumerate_cliffords(2))))
11520
```

## Common entangled states

Bell states and GHZ states have convenience constructors:

```jldoctest
julia> bell()
+ XX
+ ZZ

julia> bell(2)
+ XX__
+ ZZ__
+ __XX
+ __ZZ

julia> ghz(4)
+ XXXX
+ ZZ__
+ _ZZ_
+ __ZZ
```