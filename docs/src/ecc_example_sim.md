# [Monte Carlo simulations of noisy Clifford circuits](@id noisycircuits_mc)

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
using BenchmarkTools

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

TODO: there are Quantikz.jl bugs in the circuit visualizations

The most straightforward way to start sampling stabilizer is to set up a table of Pauli frames.

```@example 1
nframes = 200
dataqubits = code_n(code)
ancqubits = code_s(code)
regbits = ancqubits
frames = PauliFrame(nframes, dataqubits+ancqubits, regbits)
circuit = [ecirc..., scirc...]
@btime pftrajectories(frames, circuit)
```

And then you can extract the data stored in the registers:

```@example 1
frames = PauliFrame(nframes, dataqubits+ancqubits, regbits)
pftrajectories(frames, circuit)
pfmeasurements(frames)
```

To avoid runtime dispatch, it is a good practice to compactify the circuit:

```@example 1
frames = PauliFrame(nframes, dataqubits+ancqubits, regbits)
ccircuit = compactify_circuit(circuit)
@btime pftrajectories(frames, ccircuit)
```

TODO: there are spurious allocations, fix them

If you want to model Pauli errors, use:

- The helper [`PauliError`](@ref) for unbiased Pauli noise operation acting on a given qubit
- The lower level [`NoiseOp`](@ref) (for a single qubit) or [`NoiseOpAll`](@ref) (for all qubits) parameterized with a particular noise type, e.g. [`UnbiasedUncorrelatedNoise`](@ref)

TODO: example and compactification to be implemented