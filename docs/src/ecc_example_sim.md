# [ECC example with Pauli Frames](@id noisycircuits_pf_ecc_example)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using Quantikz
end
CurrentModule = QuantumClifford.Experimental.NoisyCircuits
```

Consider Steane 7-qubit code:

```@example 1
using QuantumClifford
using QuantumClifford.ECC: Steane7, naive_syndrome_circuit, encoding_circuit, parity_checks, code_s, code_n
using Quantikz

code = Steane7()
H = parity_checks(code)
```

... and the corresponding encoding circuit
```@example 1
ecirc = encoding_circuit(code)
```

... and the corresponding syndrome measurement circuit (the non-fault tolerant one)
```@example 1
scirc = naive_syndrome_circuit(code)
```

The most straightforward way to start sampling syndromes is to set up a table of Pauli frames.

```@example 1
nframes = 4
dataqubits = code_n(code)
ancqubits = code_s(code)
regbits = ancqubits
frames = PauliFrame(nframes, dataqubits+ancqubits, regbits)
circuit = [ecirc..., scirc...]
pftrajectories(frames, circuit) # run the sims
pfmeasurements(frames)          # extract the measurements
```

To avoid runtime dispatch, it is a good practice to compactify the circuit:

```@example 1
frames = PauliFrame(nframes, dataqubits+ancqubits, regbits)
ccircuit = compactify_circuit(circuit)
pftrajectories(frames, ccircuit)
pfmeasurements(frames)
```

If you want to model Pauli errors, use:

- The helper [`PauliError`](@ref) for unbiased Pauli noise operation acting on a given qubit
- The lower level [`NoiseOp`](@ref) (for a single qubit) or [`NoiseOpAll`](@ref) (for all qubits) parameterized with a particular noise type, e.g. [`UnbiasedUncorrelatedNoise`](@ref)

```@example 1
errprob = 0.1
errors = [PauliError(i,errprob) for i in 1:dataqubits]
fullcircuit = [ecirc..., errors..., scirc...]
```

And running this noisy simulation:
```@example 1
frames = PauliFrame(nframes, dataqubits+ancqubits, regbits)
pftrajectories(frames, fullcircuit)
pfmeasurements(frames)
```