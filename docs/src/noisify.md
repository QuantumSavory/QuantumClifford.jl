## Noise Insertion

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

The `noisify` API converts an ideal circuit into a noisy circuit by inserting `NoiseOp`s at configurable locations. The original circuit remains unchanged and a new noisy circuit is returned.

Using this functionality, you can apply:

- A simple noise model, where the same noise is applied to all operations.
- A structured `CircuitNoise` model, where different noise can be assigned to single-qubit gates, two-qubit gates, idle qubits, measurements, and reset operations.

### Simple Noise Model

```@example noisify
using QuantumClifford
state = one(Stabilizer, 1)

circuit = [
    sHadamard(1),
    sCNOT(1,2),
    sMZ(1,1),
    Reset(state, [2])
]

noisy = noisify(circuit, PauliNoise(1e-3, 1e-3, 1e-3))
```

This inserts noise before single-qubit gates, two-qubit gates, measurements, and reset operations.

### Structured Noise Model

```@example noisify
state = one(Stabilizer, 1)

circuit = [
    sHadamard(1),
    sCNOT(1,2),
    sMZ(1,1),
    Reset(state, [2])
]

noise_model = CircuitNoise(
    single_qubit = PauliNoise(1e-4, 1e-4, 1e-4),
    two_qubit    = PauliNoise(1e-3, 1e-3, 1e-3),
    idle_noise   = PauliNoise(1e-5, 1e-5, 1e-5),
    measurement  = PauliNoise(2e-3, 2e-3, 2e-3),
    reset        = PauliNoise(2e-3, 2e-3, 2e-3),
)

noisy_circuit = noisify(circuit, noise_model; nqubits=2)
```

To avoid double-noisification, `noisify` leaves existing `NoiseOp`s unchanged. It also does not modify `VerifyOp`s or classical operations such as `ClassicalXOR`.

### Idle Noise

When the `idle_noise` option is configured in `CircuitNoise`, noise is inserted on qubits that are not participating in any operation during a given circuit layer.

Idle noise requires the total number of qubits in the circuit to be specified.

```@example noisify
noisy_circuit = noisify(circuit, noise_model; nqubits=2)
```

If idle noise is requested and `nqubits` is not provided, an error is thrown.

### Simulation

The resulting noisy circuits can be simulated directly using the existing trajectory simulators.

```@example noisify
register = Register(one(Stabilizer, 2), falses(1))

mctrajectories(register, noisy_circuit; trajectories=100)
```
