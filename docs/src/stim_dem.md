# [Importing Stim Detector Error Models](@id stim-dem)

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

Many quantum error-correction workflows use Stim [detector error model
(`.dem`) files](https://github.com/quantumlib/Stim/blob/main/doc/file_format_dem_detector_error_model.md)
as an interchange format. A detector error model is a list of independent
classical error mechanisms: each `error(p)` instruction fires independently with
probability `p` and, when it fires, flips the *detectors* (`D#`) and *logical
observables* (`L#`) it lists.

[`read_detector_error_model`](@ref) (or [`parse_detector_error_model`](@ref) for
strings) imports such a file as a circuit -- a vector of operations -- which is
sampled by the existing [`pftrajectories`](@ref) machinery:

```julia
using QuantumClifford

circuit = read_detector_error_model("surface_code.dem")
frames = pftrajectories(circuit; trajectories=10_000)
measurements(frames)
```

## Output dimensions

For a model with ``D`` detectors and ``L`` logical observables,
`measurements(frames)` is a `trajectories × (D + L)` `Bool` matrix:

- column ``k+1`` holds detector `Dk` (columns `1:D`),
- column ``D+k+1`` holds logical observable `Lk` (columns `D+1:D+L`).

An entry is `true` when the corresponding detector or observable was flipped in
that trajectory. [`detectorview`](@ref) and [`observableview`](@ref) slice the two
blocks apart. Here is a complete deterministic example (error probabilities 0
and 1, so every trajectory is identical):

```jldoctest
julia> circuit = parse_detector_error_model("""
           error(1) D0 D2 L0
           error(0) D1
           """);

julia> circuit.n_detectors, circuit.n_logicals
(3, 1)

julia> frames = pftrajectories(circuit; trajectories=4, threads=false);

julia> measurements(frames) # trajectories × (D + L) -- here 4×4
4×4 Matrix{Bool}:
 1  0  1  1
 1  0  1  1
 1  0  1  1
 1  0  1  1

julia> collect(detectorview(circuit, frames))
4×3 Matrix{Bool}:
 1  0  1
 1  0  1
 1  0  1
 1  0  1

julia> collect(observableview(circuit, frames))
4×1 Matrix{Bool}:
 1
 1
 1
 1
```

The imported circuit consists of one no-op [`DemDeclaration`](@ref) (which makes
sure that even detectors and observables that no error ever flips are allocated
output columns) followed by one [`DetectorError`](@ref) per error mechanism:

```jldoctest
julia> parse_detector_error_model("""
           detector D0
           detector D1
           logical_observable L0
           error(0.1) D0
           error(0.2) D0 D1 L0
           """)
3-element DetectorErrorModelCircuit:
 DemDeclaration(2, 1)
 DetectorError(0.1, [1], Int64[])
 DetectorError(0.2, [1, 2], [3])
```

## Supported syntax

The parser handles the useful subset of the `.dem` grammar: `error(p)` with
`D#`/`L#` targets and `^` separators, `detector` (with optional coordinate
arguments), `logical_observable`, `shift_detectors`, nested `repeat` blocks,
instruction tags, comments, and blank lines. Instruction names are
case-insensitive. The semantics match Stim's documentation:

- each `error(p)` line is one independent Bernoulli mechanism (mechanisms with
  identical targets are *not* fused -- they remain independent, as in Stim's
  sampler);
- targets repeated within one instruction cancel out (so `error(p) D2 L0 ^ D3 L0`
  flips `D2` and `D3` but not `L0`), while still counting towards the model size;
- `shift_detectors` offsets and `repeat` blocks (which are unrolled at import
  time) only affect detector indices, never logical observable indices;
- the number of detectors is the largest mentioned or declared absolute detector
  index plus one (matching `stim.DetectorErrorModel.num_detectors`), and likewise
  for observables.

Coordinate arguments and instruction tags are parsed but discarded, and the
decomposition hints marked by `^` do not affect sampling. Anything else --
unknown instructions, malformed targets, out-of-range probabilities, unbalanced
blocks -- raises an `ArgumentError` identifying the offending line.

## API

The relevant docstrings are [`read_detector_error_model`](@ref),
[`parse_detector_error_model`](@ref), [`DetectorErrorModelCircuit`](@ref),
[`DetectorError`](@ref), [`DemDeclaration`](@ref), [`detectorview`](@ref), and
[`observableview`](@ref), all rendered in the [full API list](@ref Full-API).

```@meta
DocTestSetup = nothing
```
