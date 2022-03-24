# News

## v0.4.3

- Implement `trusted_rank` that returns `rank` for states that support it and `length` for others (e.g. `Stabilizer`).
- Implement `length` for `[Mixed]Destabilizer`.
- Clean up code repetition between `project!` and `projectX/Y/Z!`.
- More conversion constructors between different tableau types.

## v0.4.2

- `project!` does not needlessly allocate anymore on `MixedDestabilizer`. PR `#41`

## v0.4.1

- `apply_single_*` are not exported anymore, as `sX`/`sY`/`sZ` are cleaner choices.
- Faster single-qubit projections with `projectX!`, `projectY!`, `projectZ!`.
- Move circuit plotting with `Quantikz.jl` to `QuantumCliffordPlots.jl`
- Random states with zeroed phases with `phases=false`.
- Pre-compilation and inference cleanup (useful for Julia 1.8+).
- Switch from `LoopVectorization.jl` to `SIMD.jl` for fine-tuning of performance (and coincidentally, better TTFX).

## v0.4.0

- Permit whitespace separators in the `S` string macro.
- **(breaking)** `project!` now returns `anticom_index=rank` instead of `anticom_index=0` in the case of projection operator commuting with all stabilizer rows but not in the stabilizer group. If your code previously had `anticom_index!=0` checks, now you might want to use `anticom_index!=0 && anticom_index!=rank`. Conversely, treating projective measurements in general code is now much simpler.
- **(fix `#31` `b86b30e2`)** Dependent on the above, a bug fix to `Experimental.DenseMeasurement` when the measurement operator commutes with but is not in the stabilizer.
- A new `expect` function to find the expectation value of a Pauli measurement for a given stabilizer; simpler to use compared to `project!`.
- **(fix `#28` `9292333a`)** Fix a rare bug in `reset_qubits!(::MixedDestabilizer)`.

## v0.3.0

- `dot` and `logdot` for calculating dot products between `Stabilizer` states.
- Initial support for graph states, e.g., conversion to and from arbitrary `Stabilizer` state.
- **(breaking)** Implemented `Makie.jl` plotting recipes in the `QuantumCliffordPlots.jl` package, which now also stores the `Plots.jl` recipes.
- Much faster `tensor` product of states.
- **(breaking)** `CliffordColumnForm` has been completely removed. Only `CliffordOperator` now exists.
- **(breaking)** `random_singlequbitop` was removed, as it was using `CliffordColumnForm`. `random_clifford1` is a partial alternative.
- Drop `Require` to improve import times.
- Simplify internal data layout for `Stabilizer`.
- **(fix `4b536231`)** Fixed bug in `generate!` that occurs on small `IZ` Paulis.
- Some performance improvements related to allocations in `apply!`.

## v0.2.12

- `apply!` is now multi-threaded thanks to Polyester.
- Named Clifford operators with much-faster special-cased `apply!` are implemented.
- An assortment of new constructors are available to ease the conversion between data structures.
- Drop support for Julia 1.5.

## Older - unrecorded
