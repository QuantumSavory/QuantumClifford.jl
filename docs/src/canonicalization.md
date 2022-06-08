# Canonicalization operations

Different types of canonicalization operations are implemented. All of them are types of Gaussian elimination.

## [`canonicalize!`](@ref)

First do elimination on all X components and only then perform elimination on
the Z components. Based on [garcia2012efficient](@cite).
It is used in [`logdot`](@ref) for inner products of stabilizer states.

The final tableaux, if square should look like the following
![](canonicalize.png)

If the tableaux is shorter than a square, the diagonals might not reach all the way to the right.

```julia
julia> plot(canonicalize!(random_stabilizer(20,30)));
```

![](plot-canostab.png)

## [`canonicalize_rref!`](@ref)

Cycle between elimination on X and Z for each qubit. Particularly useful for
tracing out qubits. Based on [audenaert2005entanglement](@cite).
For convenience reasons, the canonicalization starts from the bottom row,
and you can specify as a second argument which columns to be canonicalized
(useful for tracing out arbitrary qubits, e.g., in [`traceout!`](@ref)).

The tableau canonicalization is done in recursive steps, each one of which results in something akin to one of these three options
![](canonicalize_rref.png)

```julia
julia> plot(canonicalize_rref!(random_stabilizer(20,30),1:30)[1]; xzcomponents=:together);
```

![](plot-rref-together.png)

## [`canonicalize_gott!`](@ref)

First do elimination on all X components and only then perform elimination on
the Z components, but without touching the qubits that were eliminated during
the X pass.
Unlike other canonicalization operations, qubit columns are reordered,
providing for a straight diagonal in each block.
Particularly useful as certain blocks of the new created matrix are
related to logical operations of the corresponding code,
e.g. computing the logical X and Z operators of a [`MixedDestabilizer`](@ref).
Based on [gottesman1997stabilizer](@cite).

A canonicalized tableau would look like the following (the right-most block does
not exist for square tableaux).
![](canonicalize_gott.png)

```julia
julia> plot(canonicalize_gott!(random_stabilizer(30))[1]; xzcomponents=:together);
```

![](plot-gottstab-together.png)

## [`canonicalize_clip!`](@ref)

Convert to the "clipped" gauge of a stabilizer state resulting in a "river" of non-identity operators around the diagonal.

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, Plots
plot(canonicalize_clip!(random_stabilizer(30)); xzcomponents=:together)
savefig("plot-clip-together.png"); nothing
```
![](plot-gottstab-together.png)

The properties of the clipped gauge are:

1. Each qubit is the left/right "endpoint" of exactly two stabilizer rows.
2. For the same qubit the two endpoints are always different Pauli operators.

This canonicalization is used to derive the [`bigram`](@ref) a stabilizer state,
which is also related to entanglement entropy in the state.

Introduced in [nahum2017quantum](@cite), with a more detailed explanation of the algorithm in Appendix A of [li2019measurement](@cite).
