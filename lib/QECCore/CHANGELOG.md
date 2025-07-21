# News

## v0.1.2 - dev

- Introduce `metacheck_matrix_x`, `metacheck_matrix_z`, and `metacheck_matrix` for CSS codes built using chain complexes and homology.
- Move the following codes from `QuantumClifford.ECC` to `QECCore`: `ReedMuller`, `RecursiveReedMuller`, `QuantumReedMuller`, `Hamming`, `Golay`, `Triangular488 `, `Triangular666 `, `Gottesman`.
- Add `Delfosse-Reichardt` codes from classical self-orthogonal `Reed-Muller` seed codes to `QECCore`.
- Add `[[4p, 2(p − 2), 4]]` Delfosse-Reichardt repetition `DelfosseReichardtRepCode` code to `QECCore`.
- Add `[[8p, 4p − 2, 3]]` Delfosse-Reichardt Generalized `[[8,2,3]]` `DelfosseReichardt823` code to `QECCore`.
- Add novel `[[n² + m²,(n - rank([C ∣ M]))² + (m − rank([C ∣ M]ᵀ))², d]]` quantum Tillich-Zémor code and `random_TillichZemor_code` code to `QECCore`.
- Introduce `QECCoreNemoExt` extension in `QECCore` for accurate `rank` computation.

## v0.1.1

- Add `AbstractDistanceAlg` for dispatch for `distance`

## v0.1.0

- First release, moving basic definitions from QuantumClifford.jl
