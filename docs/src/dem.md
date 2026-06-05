# [Importing Stim Detector Error Models](@id dem-import)

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

Many quantum error-correction workflows use [Stim](https://github.com/quantumlib/Stim)
[detector error model (`.dem`)](https://github.com/quantumlib/Stim/blob/main/doc/file_format_dem_detector_error_model.md)
files as an interchange format. A detector error model is a list of *independent
error mechanisms*; each mechanism has a probability and a set of *symptoms*
(detectors) and *frame changes* (logical observables) that it flips when it
fires.

QuantumClifford can import a `.dem` file into a circuit and sample it with the
existing fast [`pftrajectories`](@ref) Pauli-frame backend, using
[`read_detector_error_model`](@ref):

```julia
using QuantumClifford

circuit = read_detector_error_model("surface_code.dem")
frames  = pftrajectories(circuit; trajectories=10_000)
measurements(frames)
```

No new sampling engine is introduced: each `error(p) D... L...` instruction is
lowered to a [`DetectorError`](@ref) operation that, with probability `p`, flips
the listed detector and logical-observable bits of each Pauli frame.

## Output layout

[`read_detector_error_model`](@ref) returns a [`DetectorErrorModelCircuit`](@ref),
which behaves like a vector of operations but additionally records
`num_detectors` and `num_logicals`. These describe the columns of the
`trajectories × (num_detectors + num_logicals)` measurement matrix returned by
[`measurements`](@ref):

  * columns `1:num_detectors` hold the detector outcomes `D0, D1, ...`
  * columns `num_detectors .+ (1:num_logicals)` hold the observable outcomes `L0, L1, ...`

## A minimal example

Consider the small model from the file below:

```@example dem
using QuantumClifford # hide
dem = """
detector D0
detector D1
logical_observable L0
error(0.1) D0
error(0.2) D0 D1 L0
"""

# `read_detector_error_model` also accepts an `IO`, which is convenient for demos
circuit = read_detector_error_model(IOBuffer(dem))
```

It declares two detectors, one logical observable, and two independent error
mechanisms. Sampling it produces detector/logical frequencies consistent with
the declared probabilities:

```@example dem
frames = pftrajectories(circuit; trajectories=1_000_000)
m = measurements(frames)

D0 = sum(m[:, 1]) / size(m, 1)                                   # ≈ 0.1*0.9 + 0.2*0.9
D1 = sum(m[:, 2]) / size(m, 1)                                   # ≈ 0.2
L0 = sum(m[:, circuit.num_detectors + 1]) / size(m, 1)           # ≈ 0.2
(; D0, D1, L0)
```

The second mechanism `error(0.2) D0 D1 L0` always flips `D1` and `L0` together,
so those two columns are perfectly correlated.

## Supported syntax

The parser handles the commonly used subset of the `.dem` format:

  * `error(p) D... L...` error mechanisms (with `^` decomposition separators —
    only the overall set of symptoms/frame changes matters for sampling, so
    targets appearing an even number of times cancel out);
  * `detector ...` and `logical_observable ...` declarations (their coordinate
    annotations are parsed but ignored);
  * `shift_detectors ...`, including its detector-index offset, which is what
    makes loop-folded models such as Stim's repetition- and surface-code outputs
    expand correctly;
  * `repeat K { ... }` blocks, including nesting;
  * comments (`# ...`), tags (`error[tag](...)`), and blank lines.

Unsupported instructions or malformed lines raise an `ArgumentError` describing
the offending line.

## See also

  * [`read_detector_error_model`](@ref)
  * [`DetectorError`](@ref)
  * [`DetectorErrorModelCircuit`](@ref)
  * [`pftrajectories`](@ref)
