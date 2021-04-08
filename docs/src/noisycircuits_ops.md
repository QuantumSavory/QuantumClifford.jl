# [Operators in Circuit Simulations](@id noisycircuit_ops)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.Experimental.NoisyCircuits
end
CurrentModule = QuantumClifford.Experimental.NoisyCircuits
```

Available by `using QuantumClifford.Experimental.NoisyCircuits`.

## Unitary Gates

They can be specified by giving a Clifford operator tableaux and the indices on which it acts
(particularly useful for gates acting on a small part of a circuit):

```@example 1
using QuantumClifford # hide
using QuantumClifford.Experimental.NoisyCircuits # hide
SparseGate(CNOT, [2,4])
```

The Clifford operator tableaux can be completely arbitrary.
```@example 1
SparseGate(random_clifford(3), [2,4,5])
```

If the Clifford operator acts on all qubits, we do not need to specify indices.
```@example 1
DenseGate(random_clifford(5))
```

## Noisy Gates

Each gate can be followed by noise applied to the qubits on which it has acted.
This is done by wrapping the given gate into a [`NoisyGate`](@ref)

```@example 1
ε = 0.03 # X/Y/Z error probability
noise = UnbiasedUncorrelatedNoise(ε)
noisy_gate = NoisyGate(SparseGate(CNOT, [2,4]), noise)
```

In circuit diagrams the noise is not depicted, but after each application of the gate defined in `noisy_gate`, a noise operator will also be applied. The example above is of Pauli Depolarization implemented by [`UnbiasedUncorrelatedNoise`](@ref).

One can also apply only the noise operator by using [`NoiseOp`](@ref) which acts only on specified qubits. Or alternatively, one can use [`NoiseOpAll`](@ref) in order to apply noise to all qubits.

```@example 1
[NoiseOp(noise, [4,5]), NoiseOpAll(noise)]
```

## Coincidence Measurements

Global parity measurements involving single-qubit projections and classical communication are implemented with [`BellMeasurement`](@ref). One needs to specify the axes of measurement and the qubits being measured. If the parity is trivial, the circuit continues, if the parity is non-trivial, the circuit ends and reports a detected failure.
This operator is frequently used in the simulation of entanglement purification.

```@example 1
BellMeasurement([X, Y, Z], [1,3,4])
```

TODO: [`NoisyBellMeasurement`](@ref)

## Stabilizer Measurements

A measurement over one or more qubits can also be performed, e.g., a direct stabilizer measurement on multiple qubits without the use of ancillary qubits. When applied to multiple qubits, this differs from `BellMeasurement` as it performs a single projection, unlike `BellMeasurement` which performs a separate projection for every single qubit involved. This measurement is implemented in [`DenseMeasurement`](@ref) which requires a Pauli operator on which to project and the index of the classical bit in which to store the result.

TODO: example and noisy measurement, SparseMeasurement

## Reset Operations

TODO

## Gates Conditioned on Classical Bits

ΤODO: `ConditionalGate`, `DecisionGate`

TODO: Split `ConditionalGate` into quantum conditional and classical conditional