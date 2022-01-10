# News

## v0.3.0

- `dot` and `logdot` for calculating dot products between `Stabilizer` states.
- Initial support for graph states, e.g., conversion to and from arbitrary `Stabilizer` state.
- Implemented `Makie.jl` plotting recipes in the `QuantumCliffordPlots.jl` package, which now also stores the `Plots.jl` recipes.
- Much faster `tensor` product of states.
- `CliffordColumnForm` has been completely removed. Only `CliffordOperator` now exists.
- `random_singlequbitop` was removed, as it was using `CliffordColumnForm`. `random_clifford1` is a partial alternative.
- Drop `Require` to improve import times.
- Simplify internal data layout for `Stabilizer`.
- Fixed bug in `generate!` that occurs on small `IZ` Paulis.
- Some performance improvements related to allocations in `apply!`.

## v0.2.12

- `apply!` is now multi-threaded thanks to Polyester.
- Named Clifford operators with much-faster special-cased `apply!` are implemented.
- An assortment of new constructors are available to ease the conversion between data structures.
- Drop support for Julia 1.5.

## Older - unrecorded