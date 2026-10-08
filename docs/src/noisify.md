## Noise Insertion

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

The [`noisify`](@ref) API converts an ideal circuit into a noisy circuit by inserting [`NoiseOp`](@ref)'s at configurable locations. The original circuit remains unchanged and a new noisy circuit is returned.

Using this functionality, you can apply:

- A simple noise model, where the same noise is applied to all operations.
- A structured [`CircuitNoise`](@ref) model, where different noise can be assigned to single-qubit gates, two-qubit gates, idle qubits, measurements, and reset operations.

[`noisify`](@ref) does not add noise to [`VerifyOp`](@ref) or classical operations such as [`ClassicalXOR`](@ref).

### Simple Noise Model

```@example noisify
using QuantumClifford # hide

reset_state = one(Stabilizer, 1)

circuit = [
    sHadamard(1),
    sCNOT(1,2),
    sMZ(1,1),
    Reset(reset_state, [2]),
    VerifyOp(one(Stabilizer, 1), [2]),
]

noisy = noisify(circuit, PauliNoise(1e-3, 1e-3, 1e-3))
```

This applies the same noise to all supported noise locations, including single-qubit gates, two-qubit gates, measurements, reset operations, and idle qubits.


!!! note
    Idle noise is not inserted when using a simple noise model as shown above. To apply idle noise, use a [`CircuitNoise`](@ref)
    configuration with the `idle_noise` field specified.

### Structured Noise Model

A [`CircuitNoise`](@ref) object allows different noise models to be assigned to different noise operation types.

```@example noisify
noise_model = CircuitNoise(
    single_qubit = PauliNoise(1e-4, 1e-4, 1e-4),
    two_qubit    = PauliNoise(1e-3, 1e-3, 1e-3),
    idle_noise   = PauliNoise(1e-5, 1e-5, 1e-5),
    measurement  = PauliNoise(2e-3, 2e-3, 2e-3),
    reset        = PauliNoise(2e-3, 2e-3, 2e-3),
)

noisy_circuit = noisify(circuit, noise_model)
```

When the `idle_noise` field of [`CircuitNoise`](@ref) is configured, idle noise is inserted on qubits that are not involved in any operation within a circuit layer.

!!! note
    Idle noise is only applied to qubits that are part of the circuit. Qubits that exist in the underlying stabilizer state or register but never appear in any circuit operation are not considered by [`noisify`](@ref).

### Trajectory Simulation

Noisy circuits generated with [`noisify`](@ref) can be simulated using the existing trajectory-based simulators provided by QuantumClifford.

The following example demonstrates a noisy circuit simulated using Monte Carlo trajectories with [`mctrajectories`](@ref).

```@example noisify
state = one(MixedDestabilizer, 2)

v = VerifyOp(S"XX ZZ", [1,2])

mc_circuit = [
    sHadamard(1),
    sCNOT(1,2),
    v
]

noise_model = CircuitNoise(
    single_qubit = PauliNoise(1e-3, 1e-3, 1e-3),
    two_qubit = PauliNoise(1e-3, 1e-3, 1e-3),
    measurement = PauliNoise(1e-3, 1e-3, 1e-3),
)

mc_noisy_circuit = noisify(mc_circuit, noise_model)

mctrajectories(state, mc_noisy_circuit; trajectories=100)
```

Noisy circuits can also be simulated using the Pauli-frame simulator [`pftrajectories`](@ref).

```@example noisify
pf_circuit = [
    sHadamard(1),
    sCNOT(1,2),
    sMZ(1,1),
]

pf_noisy_circuit = noisify(pf_circuit, noise_model)

register = Register(state, falses(1))

pftrajectories(register, pf_noisy_circuit; trajectories=100)
```

All other simulator routines are also supported by this API