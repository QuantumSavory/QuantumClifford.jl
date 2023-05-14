# QuantumClifford.jl

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

QuantumClifford.jl is a Julia library for simulation of Clifford circuits, which are a subclass of quantum circuits that can be efficiently simulated on a classical computer.

This library uses the tableaux formalism[^1] with the destabilizer improvements[^2]. Pauli frames are supported for faster repeated simulation of noisy circuits. Various symbolic and algebraic tools for manipulating, converting, and visualizing states and circuits are also implemented. 

[^1]: [gottesman1998heisenberg](@cite)

[^2]: [aaronson2004improved](@cite)

The library consists of two main parts: Tools for working with the algebra of Stabilizer Tableaux and tools specifically for efficient Circuit Simulation.

## Stabilizer Tableau Algebra

The Stabilizer Tableau Algebra component of QuantumClifford.jl efficiently handles [pure](@ref Stabilizers) and [mixed stabilizer](@ref Mixed-Stabilizer-States) states of thousands of qubits, along with support for [sparse or dense Clifford operations](@ref Clifford-Operators) acting upon them. It provides operations such as [canonicalization](@ref Canonicalization-of-Stabilizers), [projection](@ref Projective-Measurements), [generation](@ref Generating-a-Pauli-Operator-with-Stabilizer-Generators) , and [partial traces](@ref Partial-Traces). The code is vectorized and multithreaded, offering fast, in-place, and allocation-free implementations. Tools for conversion to [graph states](@ref Graph-States) and for [visualization of tableaux](@ref Visualizations) are available.

See the [Stabilizer Tableau Algebra manual](@ref Stabilizer-Tableau-Algebra-Manual) or the curated list of [useful functions](@ref Full-API).

### Example Usage

```julia
julia> using QuantumClifford

julia> P"X" * P"Z"
-iY

julia> P"X" ⊗ P"Z"
+ XZ

julia> S"-XX
         +ZZ"
- XX
+ ZZ

julia> tCNOT * S"-XX
                 +ZZ"
- X_
+ _Z
```

## Circuit Simulation

The circuit simulation component of QuantumClifford.jl enables Monte Carlo (or symbolic) simulations of noisy Clifford circuits. It provides three main simulation methods: `mctrajectories`, `pftrajectories`, and `petrajectories`. These methods offer varying levels of efficiency, accuracy, and insight.

### Monte Carlo Simulations with Stabilizer Tableaux (`mctrajectories`)

The `mctrajectories` method runs Monte Carlo simulations using a Stabilizer tableau representation for the quantum states.

### Monte Carlo Simulations with Pauli Frames (`pftrajectories`)

The `pftrajectories` method runs Monte Carlo simulations of Pauli frames over a single reference Stabilizer tableau simulation. This approach is much more efficient but supports a smaller class of circuits.

### Symbolic Depth-First Traversal of Quantum Trajectories (`petrajectories`)

The `petrajectories` method performs a depth-first traversal of the most probable quantum trajectories, providing a fixed-order approximation of the circuit's behavior. This approach gives symbolic expressions for various figures of merit instead of just a numeric value.