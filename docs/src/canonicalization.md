# Canonicalization operations

Three different types of canonicalization operations are implemented.

All of them are types of Gaussian elimination, with different choices for how to
relate the X and Z components.

## [`canonicalize!`](@ref)

First do elimination on all X components and only then perform elimination on
the Z components. Based on [garcia2012efficient](@cite)

```julia
julia> plot(canonicalize!(random_stabilizer(20,30)));
```

![](plot-canostab.png)

## [`canonicalize_rref!`](@ref)

Cycle between elimination on X and Z for each qubit. Particularly useful for
tracing out qubits. Based on [audenaert2005entanglement](@cite).

```julia
julia> plot(canonicalize_rref!(random_stabilizer(20,30),1:30)[1]; xzcomponents=:together);
```

![](plot-rref-together.png)

## [`canonicalize_gott!`](@ref)

First do elimination on all X components and only then perform elimination on
the Z components, but without touching the qubits that were eliminated during
the X pass. Particularly useful as certain blocks of the new created matrix are
related to logical operations of the corresponding code. Based on [gottesman1997stabilizer](@cite).

```julia
julia> plot(canonicalize_gott!(random_stabilizer(20,30))[1]; xzcomponents=:together);
```

![](plot-gottstab-together.png)
