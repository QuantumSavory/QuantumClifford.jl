# QuantumClifford.jl

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

QuantumClifford.jl is a comprehensive library for the study, simulation, and manipulation of Clifford circuits and slightly non-Clifford circuits. For it, we run weekly online "office hours" [for modeling problems in quantum information science](@ref Office-Hours).

This library uses the **tableaux formalism[^1]** with the **destabilizer improvements[^2]**. Moreover, **Pauli frames** are supported for faster repeated simulation of noisy circuits. Various **symbolic and algebraic tools for manipulating, converting, and visualizing states and circuits** are also available. Finite-Clifford-rank tools for modeling circuits with **few non-Clifford gates** are provided. 

[^1]: [gottesman1998heisenberg](@cite)

[^2]: [aaronson2004improved](@cite)

The library consists of a few related components:

- Tools for working with the algebra of Stabilizer Tableaux, through data structures like [`Stabilizer`](@ref), [`MixedDestabilizer`](@ref), and [`PauliOperator`](@ref); canonicalization routines like [`canonicalize!`](@ref) and others; algebraic operations like [`project!`](@ref), partial trace, and more.

- General unitary operation acting on states: [`PauliOperator`](@ref), dense `n×n` [`CliffordOperator`](@ref), sparse small [named single-qubit and two-qubit operators](@ref Symbolic-Clifford-Operators).

- Fast circuit simulators over stabilizer states and small-Clifford-rank non-stabilizer states: Monte Carlo with [`mctrajectories`](@ref); and symbolic perturbative expansions with [`petrajectories`](@ref).

- Even faster simulator over more restricted data structures like Pauli frames with [`pftrajectories`](@ref).

- Support for small-rank non-Clifford states and arbitrary (non-unitary) Pauli channels acting on them, as well as Clifford unitaries and Pauli measurements.

- Generators of sophisticated error correcting codes, related circuits, and their evaluation.

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

The [`mctrajectories`](@ref) method runs Monte Carlo simulations using a Stabilizer tableau representation (or more general representations) for the quantum states. This simulation methods supports many operations that are not supported by the Pauli frame simulator, but it is much slower.

### Monte Carlo Simulations with Pauli Frames (`pftrajectories`)

The [`pftrajectories`](@ref) method runs Monte Carlo simulations of Pauli frames over a single reference Stabilizer tableau simulation. This approach is much more efficient but supports a smaller class of circuits.

### Symbolic Depth-First Traversal of Quantum Trajectories (`petrajectories`)

The [`petrajectories`](@ref) method performs a depth-first traversal of the most probable quantum trajectories, providing a fixed-order approximation of the circuit's behavior. This approach gives symbolic expressions for various figures of merit instead of just a numeric value.

## Non-Clifford Capabilities

The [`GeneralizedStabilizer`](@ref) type lets you run arbitrary Clifford circuits efficiently and also to run (non-unitary) Pauli Channels with a cost exponential in the number of channel applications.

## Error Correcting Codes

The `QuantumClifford.ECC` submodule can be used to generate the parity checks of many sophisticated codes, as well as to generate different styles of syndrome extraction circuits and to benchmark their performance with many different decoders.

## [Office Hours](@id Office-Hours)

Office hours are held every Friday from 12:30 – 1:30 PM Eastern Time via [Zoom](https://umass-amherst.zoom.us/j/95986275946?pwd=6h7Wbai1bXIai0XQsatNRWaVbQlTDr.1). Before joining, make sure to check the [Julia community events calendar](https://julialang.org/community/#events) to confirm whether office hours are happening, rescheduled, or canceled for the week. Feel free to bring any questions or suggestions!

## Support

QuantumClifford.jl is developed by [many volunteers](https://github.com/QuantumSavory/QuantumClifford.jl/graphs/contributors), managed at [Prof. Krastanov's lab](https://lab.krastanov.org/) at [University of Massachusetts Amherst](https://www.umass.edu/quantum/).

The development effort is supported by The [NSF Engineering and Research Center for Quantum Networks](https://cqn-erc.arizona.edu/), and
by NSF Grant 2346089 "Research Infrastructure: CIRC: New: Full-stack Codesign Tools for Quantum Hardware".

## Bounties

[We run many bug bounties and encourage submissions from novices (we are happy to help onboard you in the field).](https://github.com/QuantumSavory/.github/blob/main/BUG_BOUNTIES.md)