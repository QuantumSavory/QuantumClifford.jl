# [Importing Stim Detector Error Models](@id Importing-Detector-Error-Models)

```@meta
CurrentModule = QuantumClifford
```

[Stim](https://github.com/quantumlib/Stim) emits its noise models as detector
error model (`.dem`) files, and a lot of QEC tooling speaks that format. If you
already have such a file, [`read_detector_error_model`](@ref) turns it into a
plain QuantumClifford circuit, so you can sample it with [`pftrajectories`](@ref)
and feed the result to the decoders and analysis tools in this package — without
leaving the existing Pauli-frame backend.

```julia
using QuantumClifford

circuit = read_detector_error_model("surface_code.dem")
frames  = pftrajectories(circuit; trajectories=10_000)
measurements(frames)
```

`detector_error_model_circuit` is the same function under a second name, and if
you already have the model as a string (rather than a file) use
[`parse_detector_error_model`](@ref).

## What you get back

A `.dem` file is just a list of independent error mechanisms. Each `error(p)`
line fires with probability `p` and, when it does, flips a set of **detector**
bits (`D0, D1, …`) and **logical-observable** bits (`L0, L1, …`). The importer
turns every such line into one [`DetectorError`](@ref) operation and lets
`pftrajectories` do the sampling.

The return value, a [`DetectorErrorModelCircuit`](@ref), is a vector of those
operations (so it goes straight into `pftrajectories`) that also remembers how
many detectors and observables the model declared. After sampling,
`measurements(frames)` is a `trajectories × (num_detectors + num_observables)`
`Bool` matrix: detectors come first, observables after. As always with Pauli
frames the values are *relative* to the noiseless reference run — `true` means
"this bit was flipped in this shot".

## Worked example

Take the small model

```
detector D0
detector D1
logical_observable L0
error(0.1) D0
error(0.2) D0 D1 L0
```

Two detectors, one observable, two mechanisms:

```julia
using QuantumClifford

dem = """
detector D0
detector D1
logical_observable L0
error(0.1) D0
error(0.2) D0 D1 L0
"""

circuit = parse_detector_error_model(dem)
frames  = pftrajectories(circuit; trajectories=10^6)
m = measurements(frames)            # 10^6 × 3, columns are D0, D1, L0

vec(sum(m, dims=1)) ./ size(m, 1)   # ≈ [0.26, 0.2, 0.2]
```

`D0` sits in both mechanisms, so it flips whenever exactly one of them fires:
`0.1·(1-0.2) + 0.2·(1-0.1) = 0.26`. `D1` and `L0` only ride along with the second
mechanism, so they flip with probability `0.2` — and always together, since they
share a mechanism.

## Supported syntax

- `error(p) D… L…` — an error mechanism. The `^` separators Stim uses to hint at
  decompositions are irrelevant for sampling and are dropped; the mechanism flips
  every target listed.
- `detector D…` and `logical_observable L…` — declarations.
- `shift_detectors(coords) k` — adds `k` to the detector-index offset, so a later
  `D5` really means detector `5 + offset`. Coordinates are skipped.
- `repeat K { … }` — a block repeated `K` times, nesting allowed.
- comments (`#`), blank lines, and `[tag]`s.

A couple of things worth knowing:

- `repeat` blocks are expanded eagerly, one set of operations per iteration. This
  matches what Stim does internally, but it does mean a `repeat 100000` turns
  into a circuit with 100000 copies of the body.
- A detector or observable that is declared but never touched by an `error` still
  gets a column in the output (always `false`), so the matrix width is stable.
- Mechanisms touching more than a handful of targets fall back to a slower,
  non-compactified code path (with a warning). This essentially never happens for
  graphlike models, where mechanisms are weight two or three.

Anything the parser does not understand — an unknown instruction, an out-of-range
probability, unbalanced braces — raises an `ArgumentError` pointing at the line.

## API

See the [full API](@ref Full-API) for the docstrings of
[`read_detector_error_model`](@ref), [`detector_error_model_circuit`](@ref),
[`parse_detector_error_model`](@ref), [`DetectorErrorModelCircuit`](@ref),
[`DetectorError`](@ref), and [`DeclareMeasurementBits`](@ref).
