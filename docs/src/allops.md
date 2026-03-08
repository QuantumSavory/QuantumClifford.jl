# [Operations - Gates, Measurements, and More](@id all-operations)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using StableRNGs
    rng = StableRNG(42)
end
```

## Operations

Acting on quantum states can be performed either:

- In a "linear algebra" language where unitaries, measurements, and other operations have separate interfaces. This is an explicitly deterministic lower-level interface, which provides a great deal of control over how tableaux are manipulated. See the [Stabilizer Tableau Algebra Manual](@ref Stabilizer-Tableau-Algebra-Manual) as a primer on these approaches.
- Or in a "circuit" language, where the operators (and measurements and noise) are represented as circuit gates. This is a higher-level interface in which the outcome of an operation can be stochastic. The API for it is centered around the `apply!` function. Particularly useful for [Monte Carlo simulations](@ref noisycircuits_mc) and [Perturbative Expansion Symbolic Results](@ref noisycircuits_perturb).

In the circuit language, all operations can be applied on a state with the [`apply!`](@ref) function. Whether they are deterministic and their computational complexity is listed in the table below. A list of lower-level "linear algebra style" functions for more control over how an operation is performed is also given.

```@raw html
<style>
td > code {
    white-space: pre;
}
.content table td, .content table th {
	padding: 0 1em;
}
</style>
```

| Type | Deterministic | 𝒪(nˣ) | Low-level functions 
|:--|:-:|---|---|
|`AbstractOperation                      `|  |   |                        |
|` ├─ AbstractCliffordOperator           `|  |   |                        |
|` │   ├─ AbstractSymbolicOperator       `|  |   |                        |
|` │   │   ├─ AbstractSingleQubitOperator`|  |   |                        |
|` │   │   │   ├─ SingleQubitOperator    `|✔️ | n |                        |
|` │   │   │   ├─ sHadamard              `|✔️ | n |                        |
|` │   │   │   ├─ sId1                   `|✔️ | n |                        |
|` │   │   │   ├─ sInvPhase              `|✔️ | n |                        |
|` │   │   │   ├─ sPhase                 `|✔️ | n |                        |
|` │   │   │   ├─ sX                     `|✔️ | n |                        |
|` │   │   │   ├─ sY                     `|✔️ | n |                        |
|` │   │   │   └─ sZ                     `|✔️ | n |                        |
|` │   │   └─ AbstractTwoQubitOperator   `|  |   |                        |
|` │   │       ├─ sCNOT                  `|✔️ | n |                        |
|` │   │       ├─ sCPHASE                `|✔️ | n |                        |
|` │   │       └─ sSWAP                  `|✔️ | n |                        |
|` │   │                                 `|  |   |                        |
|` │   ├─ CliffordOperator               `|✔️ | n³|                        |
|` │   ├─ PauliOperator                  `|✔️ | n²|                        |
|` │   └─ SparseGate                     `|✔️ |kn²|                        |
|` ├─ AbstractMeasurement                `|  |   |                        |
|` │   ├─ PauliMeasurement               `|❌ | n²| [`project!`](@ref), [`projectrand!`](@ref) |
|` │   ├─ sMX                            `|❌ | n²| [`projectX!`](@ref)    |
|` │   ├─ sMY                            `|❌ | n²| [`projectY!`](@ref)    |
|` │   └─ sMZ                            `|❌ | n²| [`projectZ!`](@ref)    |
|` │                                     `|  |   |                        |
|` ├─ AbstractResetMeasurement           `|  |   |                        |
|` │   ├─ sMRX                           `|❌ | n²| [`projectX!`](@ref)   |
|` │   ├─ sMRY                           `|❌ | n²| [`projectY!`](@ref)   |
|` │   └─ sMRZ                           `|❌ | n²| [`projectZ!`](@ref)   |
|` │                                     `|  |   |                        |
|` ├─ BellMeasurement                    `|❌ | n²|                        |
|` ├─ NoiseOp                            `|❌ |  ?| [`applynoise!`](@ref)  |
|` ├─ NoiseOpAll                         `|❌ |  ?| [`applynoise!`](@ref)  |
|` ├─ NoisyGate                          `|❌ |  ?| [`applynoise!`](@ref)  |
|` ├─ Reset                              `|✔️ |kn²| [`reset_qubits!`](@ref)|
|` ├─ AbstractNonCliffordOperator        `|  |   |                        |
|` │   ├─ sT                             `|✔️ |kn²|                        |
|` │   └─ sCCZ                           `|✔️ |kn²|                        |
|` └─ AbstractPauliChannel               `|  |   |                        |
|`     ├─ PauliChannel                   `|❌ |kn²|                        |
|`     └─ UnitaryPauliChannel            `|❌ |kn²|                        |

## Details of Operations Supported by [`apply!`](@ref)

### Unitary Gates

We distinguish between symbolic gates like [`sCNOT`](@ref) that have specialized (fast) `apply!` methods (usually just for single and two qubit gates) and general tableau representation of gates like [`CliffordOperator`](@ref) that can represent any multi-qubit gate.

Predefined unitary gates are available, like [`sCNOT`](@ref), [`sHadamard`](@ref), etc.

```@example 1
using QuantumClifford # hide
using Quantikz # hide
[sCNOT(2,4),sHadamard(2),sCPHASE(1,3),sSWAP(2,4)]
```

Any arbitrary tableaux can be used as a gate too. 

They can be specified by giving a Clifford operator tableaux and the indices on which it acts
(particularly useful for gates acting on a small part of a circuit):

```@example 1
using QuantumClifford # hide
using Quantikz # hide
SparseGate(tCNOT, [2,4])
```

The Clifford operator tableaux can be completely arbitrary.
```@example 1
SparseGate(random_clifford(3), [2,4,5])
```

If the Clifford operator acts on all qubits, we do not need to specify indices, just use the operator.

### Noisy Gates

Each gate can be followed by noise applied to the qubits on which it has acted.
This is done by wrapping the given gate into a [`NoisyGate`](@ref)

```@example 1
ε = 0.03 # X/Y/Z error probability
noise = UnbiasedUncorrelatedNoise(ε)
noisy_gate = NoisyGate(SparseGate(tCNOT, [2,4]), noise)
```

In circuit diagrams the noise is not depicted, but after each application of the gate defined in `noisy_gate`, a noise operator will also be applied. The example above is of Pauli Depolarization implemented by [`UnbiasedUncorrelatedNoise`](@ref).

One can also apply only the noise operator by using [`NoiseOp`](@ref) which acts only on specified qubits. Or alternatively, one can use [`NoiseOpAll`](@ref) in order to apply noise to all qubits.

```@example 1
[NoiseOp(noise, [4,5]), NoiseOpAll(noise)]
```

The machinery behind noise processes and different types of noise is detailed in [the section on noise](@ref noise)

### Coincidence Measurements

Global parity measurements involving single-qubit projections and classical communication are implemented with [`BellMeasurement`](@ref). One needs to specify the axes of measurement and the qubits being measured. If the parity is trivial, the circuit continues, if the parity is non-trivial, the circuit ends and reports a detected failure.
This operator is frequently used in the simulation of entanglement purification.

```@example 1
BellMeasurement([sMX(1), sMY(3), sMZ(4)])
```

There is also [`NoisyBellMeasurement`](@ref) that takes the bit-flip probability of a single-qubit measurement as a third argument.

### Stabilizer Measurements

A measurement over one or more qubits can also be performed, e.g., a direct stabilizer measurement on multiple qubits without the use of ancillary qubits. When applied to multiple qubits, this differs from `BellMeasurement` as it performs a single projection, unlike `BellMeasurement` which performs a separate projection for every single qubit involved. This measurement is implemented in [`PauliMeasurement`](@ref) which requires a Pauli operator on which to project and the index of the classical bit in which to store the result. Alternatively, there are [`sMX`](@ref), [`sMZ`](@ref), [`sMY`](@ref) if you are measuring a single qubit.

```@example 1
[PauliMeasurement(P"XYZ", 1), sMZ(2, 2)]
```

### Reset Operations

The [`Reset`](@ref) operations lets you trace out the specified qubits and set their state to a specific tableau.

```@example 1
new_state = random_stabilizer(3)
qubit_indices = [1,2,3]
Reset(new_state, qubit_indices)
```

It can be done anywhere in a circuit, not just at the beginning.

### Non-Clifford Gates

The gates [`sT`](@ref) (T gate, π/8 phase rotation) and [`sCCZ`](@ref) (controlled-controlled-Z) are non-Clifford gates. They work with [`PureGeneralizedStabilizer`](@ref) states and are simulated via [`emtrajectories`](@ref), which performs end-of-circuit Z-basis measurements using the sum-over-Cliffords decomposition.

```@example 1
circuit = [sHadamard(1), sT(1), sHadamard(1)]
result = emtrajectories(circuit; trajectories=100, delta=0.1)
meas = measurements(result)
size(meas) # 100 samples × 1 qubit
```

### Pauli Channels

[`UnitaryPauliChannel`](@ref) and [`PauliChannel`](@ref) represent non-Clifford channels as weighted sums of Pauli operators. They work with [`GeneralizedStabilizer`](@ref) states via [`apply!`](@ref) and [`mctrajectories`](@ref).

Predefined channels include [`pcT`](@ref) (the T gate as a Pauli channel), [`pcPhase`](@ref)`(ϕ)` (arbitrary phase rotation), and [`pcRx`](@ref)`(θ)` (X rotation).

```@example 1
using QuantumClifford.PauliChannelNonClifford: pcT # hide
state = GeneralizedStabilizer(S"-X")
apply!(state, pcT)
```