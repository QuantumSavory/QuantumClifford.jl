# News

## v0.6.2 - 2022-09-22

- Performance improvements: fine tuning of `:constprop`, eliminating a handful of common causes of dynamic dispatch.
- Drop support for Julia 1.6.

## v0.6.1 - 2022-09-05

- **(fix)** Fix a bug in the unexported `projectremoverand!` that occurred due to the introduction of `Tableau`.

## v0.6.0 - 2022-08-22

- **(breaking)** Split the `Stabilizer` object into a `Stabilizer` that semantically represents a state, and a general `Tableau` that does not carry such an interpretation. `Stabilizer` uses `Tableau` internally. `stab.xzs` and `stab.phases` property access would now fail. Use `tab(stab)` to get the tableau object and `phases(stab)` to get the phases.
- Simplify type parameters.

## v0.5.9 - 2022-08-17

- `Register` now works with `traceout!` and `project*rand!`.

## v0.5.8 - 2022-08-15

- Implement `projectremoverand!` which besides performing a projective measurement like `projectrand!` also removes the measured qubit from the tableau, returning a smaller tableau. Not yet exported in public API.

## v0.5.7 - 2022-07-24

- **(fix)** `apply!(S"XXX", P"X", [1])` and similar sparse Pauli applies were giving wrong results.
- Significant speedup of `petrajectories` thanks to an order of magnitude speedup in `applynoise_branches(...,::UnbiasedUncorrelatedNoise)`.
- Expanding test suite, including base functions, `Experimental.NoisyCircuits`, and others. Re-establishing tests of alternative big-packing.

## v0.5.6 - 2022-07-20

- **(fix)** Bug-fixes to `PauliMeasurement`, `Reset`, and `Register`, detected by JET.

## v0.5.5 - 2022-07-05

- **(breaking fix)** `CliffordOperator` constructor called on a square tableau occasionally resulted in incorrect results due to row-reordering during cannonicalization.
- Continuing static analysis fixes thanks to JET.
- Optimization of `canonicalize_clip!`, now resulting in much fewer allocations.

## v0.5.4 - 2022-07-03

- Start using `JET.jl` for static analysis during CI.
- The `MixedDestabilizer` constructor now accepts over redundant tableaux (tableaux with redundant rows).
- Resolved multiple method ambiguities and started testing for them with `Aqua.jl` in CI.

## v0.5.3 - 2022-06-11

- **(fix `#60`)** `enumerate_clifford` was broken due to an earlier refactor in `rowswap!`

## v0.5.2 - 2022-06-08

- Implement the clipped gauge canonicalization `canonicalize_clip!` and related functions.
- Implement `entanglement_entropy`.

## v0.5.1 - 2022-06-02

- **(fix `#57` `c83f85f`)** - Graph vertex indices and qubit indices in tableaux were mismatched.

## v0.5.0 - 2022-05-19

- **(breaking)** Rename all pre-defined tableaux to have a `t` prefix. e.g., `CNOT`→`tCNOT`, in order to distinguish them from "symbolic" operators like `sCNOT`.
- **(breaking)** Rename `CliffordId` to `tId1` to match the naming style of `sId1`.
- Implement `enumerate_cliffords`, `enumerate_phases`, `symplecticGS` used for the enumeration of all Clifford operations over given number of qubits.
- Implement convertors from symbolic operators to dense tableau operators: `CliffordOperator(sCNOT)` now gives `tCNOT` which acts equivalently to `sCNOT(1,2)`.
- Implement `project[X|Y|Z]rand!` as a simpler interface to `project!` with automatic randomization of measurement phases.
- Implement `sMX`/`sMY`/`sMZ` symbolic measurement operations that can be used with `apply!`. Use `projectrand!` internally.
- **(breaking)** Cleanup of `NoisyCircuits`
  - The experimental module `NoisyCircuits` now supports only `MixedDestabilizer` and `Register`.
  - `Register` is moved out of `NoisyCircuits`. Used with `sMX`/etc to store classical bit results during circuit evolution.
  - `SparseGate` is moved out of `NoisyCircuits`.
  - `Reset` is moved out of `NoisyCircuits`.
  - `DenseMeasurement` is renamed `PauliMeasurement` and moved out of `NoisyCircuits`.
  - `DenseGate` is removed (just use any dense CliffordOperator).
  - `SparseMeasurement` is removed (just use `sMX`, `sMY`, `sMZ`). Due to this we lost the functionality of measuring more than one but less than all qubits.
  - `applyop!` is renamed to `applywstatus!` and simplified.
  - `applyop_branches` is renamed to `applybranches` and simplified.

## v0.4.3 - 2022-03-24

- Implement `trusted_rank` that returns `rank` for states that support it and `length` for others (e.g. `Stabilizer`).
- Implement `length` for `[Mixed]Destabilizer`.
- Clean up code repetition between `project!` and `projectX/Y/Z!`. Issue `#40`
- More conversion constructors between different tableau types.
- Implement pre-compilation examples (again) for julia 1.9+.
- `generate!` does not needlessly allocate anymore, helping `project!` on `Stabilizer`. Issue `#39`

## v0.4.2 - 2022-03-22

- `project!` does not needlessly allocate anymore on `MixedDestabilizer`. Issue `#39`, PR `#41`

## v0.4.1 - 2022-03-21

- `apply_single_*` are not exported anymore, as `sX`/`sY`/`sZ` are cleaner choices.
- Faster single-qubit projections with `projectX!`, `projectY!`, `projectZ!`.
- Move circuit plotting with `Quantikz.jl` to `QuantumCliffordPlots.jl`
- Random states with zeroed phases with `phases=false`.
- Pre-compilation and inference cleanup (useful for Julia 1.8+).
- Switch from `LoopVectorization.jl` to `SIMD.jl` for fine-tuning of performance (and coincidentally, better TTFX).

## v0.4.0 - 2022-01-27

- Permit whitespace separators in the `S` string macro.
- **(breaking)** `project!` now returns `anticom_index=rank` instead of `anticom_index=0` in the case of projection operator commuting with all stabilizer rows but not in the stabilizer group. If your code previously had `anticom_index!=0` checks, now you might want to use `anticom_index!=0 && anticom_index!=rank`. Conversely, treating projective measurements in general code is now much simpler.
- **(fix `#31` `b86b30e2`)** Dependent on the above, a bug fix to `Experimental.DenseMeasurement` when the measurement operator commutes with but is not in the stabilizer.
- A new `expect` function to find the expectation value of a Pauli measurement for a given stabilizer; simpler to use compared to `project!`.
- **(fix `#28` `9292333a`)** Fix a rare bug in `reset_qubits!(::MixedDestabilizer)`.

## v0.3.0 - 2022-01-10

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

## v0.2.12 - 2021-12-32

- `apply!` is now multi-threaded thanks to Polyester.
- Named Clifford operators with much-faster special-cased `apply!` are implemented.
- An assortment of new constructors are available to ease the conversion between data structures.
- Drop support for Julia 1.5.

## Older - before 2021-10-28 unrecorded
