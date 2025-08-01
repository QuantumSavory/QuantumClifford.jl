# Planned changes for v1.0.0:

- `permute` will be a wrapper around to `QuantumInterface.permutesystems`. Documentation for `permute!` would be similarly updated
- reworking the rest of `NoisyCircuits` and moving it out of `Experimental`

# News

## v0.10.1-dev

- Add novel `[[n² + m²,(n - rank([C ∣ M]))² + (m − rank([C ∣ M]ᵀ))², d]]` quantum Tillich-Zémor `random_TillichZemor_code` codes to `QECCore` and introduce `QECCoreNemoExt` for accurate matrix `rank` computation.
- The D-dimensional Surface and Toric codes are now implemented using `Oscar`'s chain complexes and `GF2` homology in the ECC submodule.
- Add `[[2n², 2k², d]]` and `[[(n - k)² + n², k², d]]` La-cross codes via `Hecke`'s univariate polynomial ring in the ECC submodule.
- Implementing `apply_inv!` for direct application of the inverse of a given gate.
- The lifted product code constructor `LPCode` now supports non-commutative group algebras by appropriate switching left/right representations — particularly useful now that there is also an `Oscar` extension, which provides many non-abelian group constructors.
- Introduce `metacheck_matrix_x`, `metacheck_matrix_z`, and `metacheck_matrix` for CSS codes built using chain complexes and homology.
- `ReedMuller`, `RecursiveReedMuller`, and `QuantumReedMuller` are moved to `QECCore` from `QuantumClifford.ECC`.
- The phase storage type can now be parameterized, instead of hardcoded to UInt8.
- Drop support for Julia 1.9.
- Add convenience wrappers for code families: `honeycomb_color_codes` and `haah_cubic_codes` via `LPCode` construction.
- Some codes are moved to `QECCore` from `QuantumClifford`, including `Hamming`, `Golay`, `Triangular488 `, `Triangular666 `, `Gottesman`.
- Add `Delfosse-Reichardt` codes from classical self-orthogonal `Reed-Muller` seed codes to `QECCore`.
- Add `[[4p, 2(p − 2), 4]]` Delfosse-Reichardt repetition `DelfosseReichardtRepCode` code to `QECCore`.
- Add `[[8p, 4p − 2, 3]]` Delfosse-Reichardt Generalized `[[8,2,3]]` `DelfosseReichardt823` code to `QECCore`.
- Add Quantum Tanner graph product codes: general and cyclic `Q(G₁×G₂)` codes (Tanner graphs `G₁`, `G₂`) to `QECCore`.
- Add classical Gallager's LDPC code to `QECCore`.

### Private API

These changes affect internal implementation details - external packages should not rely on these!
- Add `check_repr_regular_linear` to verify `F`-linear regular representation for group algebra elements.
- Add `[2ʳ-1, 2ʳ-1-r, 3]` Hamming code.

## v0.10.0 - 2025-07-02

- **(fix)** The gates `SQRTY`, `CXYZ`, `CZYX` were computing phases incorrectly when acting on `I` stabilizers.
- **(fix)** Paulis with imaginary phases had their phases incorrectly tracked.
- **(fix)** `rowdecompose` was not accounting for the phase of the input Pauli string, leading to potential errors in non-Clifford functionality.
- Various stabilizer data structures can now be constructed out of sequences of Paulis.
- Most canonicalization routines can now act on sequences of Paulis.
- `project*!` family of functions now work on more tableau types.
- `QuantumCliffordJuMPExt` is a new extension depending on the JuMP solver suite. It enables computing minimum distance for quantum LDPC codes via Mixed Integer Programming (MIP) using `JuMP`.
- `QuantumCliffordOscarExt` is a new extension depending on the Oscar CAS. It provides more convenient ways to construct two-block-group-algebra codes and more.
- Improvements to `GeneralizedStabilizer` API including support for products and tensor product for many Clifford and non-Clifford operators and states.
- `expect` is now implemented for `GeneralizedStabilizer`.
- `projectrand!` is now implemented for `GeneralizedStabilizer`.
- `QuantumClifford` now depends on `QECCore`. `QECCore` is a new package separate from `QuantumClifford`, which specifies the interfaces error correction codes, together with declaring a number of the more essential ones. Naming changes:
  - **(breaking)** `parity_checks_x` -> `parity_matrix_x` and `parity_checks_z` -> `parity_matrix_z`
  - `AbstractECC` -> `AbstractQECC`
  - `ClassicalCode` -> `AbstractCECC`
  - Some codes are moved to `QECCore` from `QuantumClifford`, including `Toric`, `Surface`, `RepCode`, `CSS`, `Shor9`, `Steane7`, `Cleve8`, `Perfect5`, `Bitflip3`.
- **(breaking)** `StabMixture` was renamed to `GeneralizedStabilizer`.
- Constructing a `Destabilizer` out of a full-rank `Stabilizer` does not require a canonicalization anymore, i.e. `stabilizerview(Destabilizer(s))==s` is guaranteed.
- The `maximally_mixed` function is now available for creating maximally mixed multi-qubit states.
- Much faster indexing and slicing of `PauliOperator`.

## v0.9.19 - 2025-04-08

- `permutesystems` and `permutesystems!` are no implemented, deprecating `permute` and `permute!`

## v0.9.18 - 2025-02-19

- Fixes for rare crashes in the python BP decoders.
- Less repetitive error messages.

## v0.9.17 - 2025-02-18

- New memory structure and matrix inversion function for `random_destabilizer`, to reduce allocations and speed up repeated generation of many random destabilizers.
- Improvements to allocations in `apply!`

## v0.9.16 - 2024-12-29

- 100× faster unbiased `random_pauli`.
- Enhancements to `GF(2)` Linear Algebra: unexported, experimental `gf2_row_echelon_with_pivots!`, `gf2_nullspace`, `gf2_rowspace_basis`.

## v0.9.15 - 2024-12-22

- `pftrajectories` now supports fast multiqubit measurements with `PauliMeasurement` in addition to the already supported single qubit measurements `sMX/Z/Y` and workarounds like `naive_syndrome_circuit`.

## v0.9.14 - 2024-11-03

- **(fix)** `affectedqubits()` on `sMX`, `sMY`, and `sMR*`
- **(fix)** restrictive type-assert in `MixedDestabilizer` failing on views of tableaux
- Implementing additional named two-qubit gates: `sSQRTXX, sInvSQRTXX, sSQRTYY, sInvSQRTYY`

## v0.9.13 - 2024-10-30

- New error-correction group theory tools:
    - `canonicalize_noncomm` function to find a generating set with minimal anticommutivity
    - `SubsystemCodeTableau` data structure to represent the output of `canonicalize_noncomm`
    - `commutify` function to find a commutative version of a non-commutative set of Paulis with minimal changes
    - `matroid_parent` to, for set of Paulis that doesn't represent a state, find a version
    that does.
- Implementing additional named two-qubit gates: `sSWAPCX, sInvSWAPCX, sCZSWAP, sCXSWAP, sISWAP, sInvISWAP, sSQRTZZ, sInvSQRTZZ`

## v0.9.12 - 2024-10-18

- Minor compat fixes for julia 1.11 in the handling of `hgp`

## v0.9.11 - 2024-09-27

- `hcat` of Tableaux objects
- `QuantumReedMuller` codes added to the ECC module
- **(breaking)** change the convention for how to provide a representation function in the constructor of `LPCode` -- strictly speaking a breaking change, but this is not an API that is publicly used in practice

## v0.9.10 - 2024-09-26

- The lifted product class of quantum LDPC codes is implemented in the ECC submodule.
- **(fix)** `ECC.code_s` now gives the number of parity checks with redundancy. If you want the number of linearly independent parity checks, you can use `LinearAlgebra.rank`.
- Implementing many more named single-qubit gates following naming convention similar to the stim package in python.
- **(fix)** Bug fix to the `parity_checks(ReedMuller(r, m))` of classical Reed-Muller code (it was returning generator matrix).
- `RecursiveReedMuller` code implementation as an alternative implementation of `ReedMuller`.


## v0.9.9 - 2024-08-05

- `inv` is implemented for all Clifford operator types (symbolic, dense, sparse).

## v0.9.8 - 2024-08-03

- New group-theoretical tools:
    - `groupify` to generate full stabilizer group from generating set
    - `minimal_generating_set` function to find the minimal generating set of a set of operators
    - `pauligroup` to generate the full Pauli group of a certain number of qubits
    - `normalizer` to generate all Paulis that commute with a set of Paulis
    - `centralizer` to find a subset of a set of Paulis such that each element in the subset commutes with each element in the set
    - `contractor` to find a subset of Paulis in a tableau that have an identity operator on a certain qubit
    - `delete_columns` to remove the operators corresponding to a certain qubit from all Paulis in a Stabilizer
- `PauliError` can now encode biased noise during Pauli frame simulation, i.e. one can simulate only X errors, or only Y errors, or only Z errors, or some weighted combination of these.

## v0.9.7 - 2024-07-23

- **(fix `#320`)** Fix a serious correctness bug in the SIMD implementation of Pauli string multiplication (affects the correctness of canonicalization and traceout for tableaux bigger than ~500 qubits; does not affect symbolic gates or Pauli frame simulations of any scale)

## v0.9.6 - 2024-07-12

- `inv` implementation for single-qubit "symbolic" Clifford operators (subtypes of `AbstractSingleQubitOperator`).

## v0.9.5 - 2024-07-04

- Implementation of random all-to-all and brickwork Clifford circuits and corresponding ECC codes.

## v0.9.4 - 2024-06-28

- Addition of a constructor for concatenated quantum codes `Concat`.
- Addition of multiple unexported classical code constructors.
- Gate errors are now conveniently supported by the various ECC benchmark setups in the `ECC` module.
- Significant improvements to the low-level circuit compiler (the sumtype compactifier), leading to faster Pauli frame simulation of noisy circuits.
- Bump `QuantumOpticsBase.jl` package extension compat bound.
- **(fix)** Remove printing of spurious debug info from the PyBP decoder.
- **(fix)** Failed compactification of gates now only raises a warning instead of throwing an error. Defaults to slower non-compactified gates.

## v0.9.3 - 2024-04-10

- **(fix)** One of `random_pauli`'s methods was disregarding the error probability and had incorrect kwarg defaults.

## v0.9.2 - 2024-04-08

- The ECC module now has access to an iterative bitflip decoder thanks to `LDPCDecoders.jl`.
- Provide more configuration options in the `PyBeliefProp` decoders.
- **(fix)** The belief prop decoder from LDPCDecoders was counting iterations incorrectly.

## v0.9.1 - 2024-03-31

- Implemented `iscss` function to identify whether a given code is known to be a CSS (Calderbank-Shor-Steane) code.
- Added the classical Reed-Muller code in the ECC module.
- Added the surface code to the ECC module.

## v0.9.0 - 2024-03-19

- **(breaking)** The defaults in `random_pauli` are now `realphase=true` and `nophase=true`.
- **(breaking)** The convention for for flip probability in `random_pauli`.
- **(breaking)** The convention for noise probability in `UnbiasedUncorrelatedNoise` changed. The input number is the total probability for an error to occur.
- Implement an inplace `random_pauli!`, a non-allocating alternative to `random_pauli`.
- Significant improvement in the performance of the ECC decoder pipeline (but many low-hanging fruits still remain).

## v0.8.21 - 2024-03-17

- Implemented the Gottesman code family, also known as [[2^j, 2^j - j - 2, 3]] quantum Hamming codes.
- Bump the `PyQDecoders` dependency, switching to `PythonCall` behind the scenes for reliability.
- Bump the `LDPCDecoders` dependency.

## v0.8.20 - 2024-01-22

- Significant additions to the `ECC` submodule, with constructors for a few new codes (`Toric` and  generic `CSS`); incorporating many syndrome decoding algorithms (thanks to the `PyQDecoders.jl` and `LDPCDecoders.jl` packages); and providing a convenient API for evaluating code performance in different settings through the new `evaluate_decoder` function.

## v0.8.19 - 2023-12-16

- Bumping up the lower bounds of many dependencies and adding lower-bound compatibility checks to CI.

## v0.8.18 - 2023-11-22

- `ECC.faults_matrix` detects and warns when encountery codes with redundant checks.

## v0.8.17 - 2023-10-17

- **(fix)** Some `affectedqubits` methods were returning single integers instead of a one-tuple.
- The non-public `ECC` module has seen a few improvements (a `naive_encoding_circuit` implementation and a few new codes), as well as some breaking changes to API.

## v0.8.16 - 2023-08-31

- **(breaking)** Changes to the circuit generation in the non-public ECC module. Adding a Shor syndrome measurement circuit.
- Added a `ClassicalXOR` gate, for now supporting only `PauliFrame`s.

## v0.8.15 - 2023-08-16

- Initial support for GPU accelerated circuit simulation (with CUDA).
- Minor documentation fixes around `phases` and a workaround for Makie plotting regression.

## v0.8.14 - 2023-07-19

- Circuit plotting with Quantikz from inside other modules (common with Pluto) showed wrong names for gates due to how we were serializing the names. It is now fixed.

## v0.8.13 - 2023-07-18

- **(fix)** There was a bug with incorrect scaling for `Operator(::CliffordOperator)` conversions.
- A few more features to the `ECC` module's circuit generation routines.
- Quantikz circuit plotting improvements to `CliffordOperator` and `s*CY` and `sYC*`.

## v0.8.12 - 2023-07-12

- Initial implementation of non-Clifford simulation (mainly for circuits that are slightly non-Clifford, e.g. containing T gates). See `StabMixture`, `PauliChannel`, `UnitaryPauliChannel`, and `pcT`.
- `embed` implemented for `PauliOperator` and `PauliChannel`.
- Various convenience constructors that transform a tableaux or an operator into a `Ket` or `Operator` from `QuantumOptics.jl`. Use the constructors directly like `Ket(::Stabilizer)`, `Operator(::CliffordOperator)`, etc.

## v0.8.11 - 2023-07-10

- `petrajectories`, for (potentially symbolic) perturbative expansions of the result of a circuit, is moved out of `Experimental` into the public part of the interface. The underlying `petrajectory` is not made public yet due to the ad-hoc low-level return conventions for it.
- `mctrajectory` and `petrajectory` can now optionally report the end state of each trajectory, not just the circuit status (i.e. "success", "detected failure", etc).
- Internally we now use a trait system to distinguish deterministic from non-deterministic operations. Not part of the public API yet.

## v0.8.10 - 2023-07-05

- Remove Polyester.jl multithreading, leading to simpler and better compiled single-thread code. Now single-thread performance is much higher. Multithreading was generally not useful at the inner loops where it was deployed.
- Implement `fastcolumn` and `fastrow`, which transform a tableau into a memory layout optimized for column or row operations.

## v0.8.9 - 2023-07-04

- In the unexported experimental ECC module:
    - we now implement `fault_matrix` which gives the mapping between single-qubit physical errors and logical observables.
    - `MixedDestabilizer` and `Stabilizer` now have constructors when acting on an ECC object.
- `stab_to_gf2` now works with Pauli operators as well.

## v0.8.8 - 2023-06-23

- Bump `QuantumInterface` compat.

## v0.8.7 - 2023-06-22

- Better UX and threading support for `pftrajectories`.
- `affectedqubits` now more consistently returns tuples instead of vectors.
- Private `affectedbits` is now implemented.
- Many operation constructors now throw an error if they are given negative qubit indices.

## v0.8.6 - 2023-06-20

- Fixes to Quantikz circuit plotting of empty circuits and `PauliOperator`

## v0.8.5 - 2023-06-13

- Internal helper method `mul_right!` is now available for right Pauli inplace multiplication.
- Implemented `sMRZ` to reset single qubits to the |0⟩ (and respectively `sMRX` and `sMRY`).

## v0.8.4 - 2023-06-11

- Bump `QuantumInterface` to 0.2.0.

## v0.8.3 - 2023-06-10

- Improvements to printing and documentation.

## v0.8.2 - 2023-05-22

- Initial, experimental, unexported helper functions for work with error correcting codes. Might change at any time.

## v0.8.1 - 2023-05-17

- **(fix)** Fix to Quantikz circuit plotting functions.

## v0.8.0 - 2023-05-16

- **(breaking)** Deprecated the `QuantumCliffordPlots` library and moved that functionality to 3 package extensions, one for each of `Quantikz`, `Plots`, and `Makie` visualizations.
- **(breaking)** Set minimum requirements to Julia 1.9
- Initial implementation of Pauli frame simulations (with `pftrajectories` and `PauliFrame`)
- Initial support for "sum type" gates for much faster dispatch (with `compactify_circuit`)
- **(fix)** Fixes to print/show formatting.

## v0.7.2 - 2023-04-09

- Non-standard string literals `P`, `T`, `S`, and `C` (e.g. `P"X"` used to create a Pauli operator) are now not cached at compile time. Before this version `f() = P"X"` would have resulted in `f() === f()`, while now that statement would return false. Change made given that the objects created by these literals are mutable which can lead to bugs where a local variable seems to "remember" how it is being modified at each execution of a function.

## v0.7.1 - 2023-04-08

- Better printing of tableaux and operators.

## v0.7.0 - 2023-03-31

- **(breaking)** Switched the convention for results from `apply!(Register(...), sMX)` (and corresponding for `sMY` and `sMZ`). Previously the +1 eigenvector and -1 eigenvector were mapped to 1 and 0 as measurement results. This is inconvenient for boolean math and contrary to typical convection. Moving to standard convention now (+1 eigenvector corresponds to result 0). Most convenient because `eigenvalue = (-1)^measurement result`.
- The clipped gauge (`canonicallize_clip!`) now works on mixed (i.e. incomplete, i.e. non-square) tableaux

## v0.6.7 - 2022-12-27

- **(fix)** Fixed bug in `CliffordOperator(AbstractTwoQubitOperator)`.
- Fixes to inference failures (detected by JET).
- Significant test coverage increase.
- Significant improvement in TTFX thanks to fixes to Polyester.jl (for Julia 1.9).
- Stabilizing a few features, moving out of `Experimental.NoisyCircuits`
    - `applynoise!`
    - `affectedqubits`
    - `NoisyGate`
    - `NoiseOp` and `NoiseOpAll`
    - `UnbiasedUncorrelatedNoise`

## v0.6.6 - 2022-12-24

- `random_destabilizer(rank,nb_of_qubits)` now exists to provide a random `MixedDestabilizer` instance of a given `rank`.
- Stabilizing a few features, moving out of `Experimental.NoisyCircuits`
    - `BellMeasurement`
    - `mctrajectories`, `mctrajectory!`, `applywstatus!`
    - `CircuitStatus`
- `colpermute!` was turned into `Base.permute!`. It was not documented or used previously.
- `check_allrowscommute` is not exported by default anymore. It was not documented or used previously.

## v0.6.5 - 2022-12-23

- Minor API adjustments in preparation for releasing `BPGates.jl`.

## v0.6.4 - 2022-12-11

- Move the API declarations to `QuantumInterface.jl`, a mostly namespacing package shared with `QuantumOptics.jl`.

## v0.6.3 - 2022-09-24

- Bring the v0.6.2 performance improvements to Julia 1.7 and 1.6. Switch from using `constprop` to explicit static dispatch with `@valbooldispatch`.
- Support Julia 1.6 again.

## v0.6.2 - 2022-09-22

- Performance improvements: fine tuning of `constprop`, eliminating a handful of common causes of dynamic dispatch.
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
- Expanding test suite, including base functions, `Experimental.NoisyCircuits`, and others. Re-establishing tests of alternative bit-packing.

## v0.5.6 - 2022-07-20

- **(fix)** Bug-fixes to `PauliMeasurement`, `Reset`, and `Register`, detected by JET.

## v0.5.5 - 2022-07-05

- **(breaking fix)** `CliffordOperator` constructor called on a square tableau occasionally resulted in incorrect results due to row-reordering during canonicalization.
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
