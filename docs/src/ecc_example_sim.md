# [ECC example with Pauli Frames](@id noisycircuits_pf_ecc_example)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using Quantikz
end
```

!!! warning "The documentation is incomplete"
    Waiting for a better documentation than the small example below.
    Check out also the page on [ECC performance evaluators](@ref ecc_evaluating)


Consider Steane 7-qubit code:

```@example 1
using QuantumClifford
using QuantumClifford.ECC: Steane7, naive_syndrome_circuit, naive_encoding_circuit, parity_checks, code_s, code_n
using Quantikz

code = Steane7()
H = parity_checks(code)
```

... and the corresponding encoding circuit
```@example 1
ecirc = naive_encoding_circuit(code)
```

... and the corresponding syndrome measurement circuit (the non-fault tolerant one)
```@example 1
scirc, _ = naive_syndrome_circuit(code)
```

The most straightforward way to start sampling syndromes is to set up a table of Pauli frames.

```@example 1
circuit = [ecirc..., scirc...]
nframes = 4
frames = pftrajectories(circuit; trajectories=nframes) # run the sims
pfmeasurements(frames)                                 # extract the measurements
```

The [`pftrajectories`](@ref) function is multithreaded.
If you want more low-level control over these Pauli frame simulations, check out the [`PauliFrame`](@ref) structure,
the other methods of [`pftrajectories`](@ref), and the circuit compactifaction function [`compactify_circuit`](@ref).

If you want to model Pauli errors, use:

- The helper [`PauliError`](@ref) for unbiased Pauli noise operation acting on a given qubit
- The lower level [`NoiseOp`](@ref) (for a single qubit) or [`NoiseOpAll`](@ref) (for all qubits) parameterized with a particular noise type, e.g. [`UnbiasedUncorrelatedNoise`](@ref)

```@example 1
errprob = 0.1
errors = [PauliError(i,errprob) for i in 1:code_n(code)]
fullcircuit = [ecirc..., errors..., scirc...]
```

And running this noisy simulation:
```@example 1
frames = pftrajectories(fullcircuit; trajectories=nframes)
pfmeasurements(frames)
```
## Exporting a Stim detector error model

For interoperability with [Stim](https://github.com/quantumlib/Stim)-based workflows and
external decoders, a code-capacity detector error model can be generated from any code
that provides [`parity_checks`](@ref) and [`faults_matrix`](@ref), and written out in
Stim's `.dem` text format. Each enabled single-qubit Pauli fault (with probabilities
given independently by `px`, `py`, `pz`) becomes one `error(p)` instruction, whose `D`
targets are the parity checks flipped by that fault and whose `L` targets are the logical
observables flipped by it. Detector `Dᵢ` corresponds to row `i+1` of `parity_checks(code)`
and observable `Lⱼ` to row `j+1` of `faults_matrix(code)`. The output is deterministic,
so the files are easy to test and diff.

```@example 1
using QuantumClifford.ECC # hide
dem = detector_error_model(Steane7(); px=1e-3, py=0.0, pz=1e-3)
```

```@example 1
write_detector_error_model(stdout, dem)
```

The same model can be written directly to a file with
`write_detector_error_model("steane7.dem", dem)`, which Stim can then parse:

```python
>>> import stim
>>> dem = stim.DetectorErrorModel.from_file("steane7.dem")
>>> dem.num_detectors, dem.num_observables
(6, 2)
```
