# News

## v0.1.2 - 2025-12-31

- Introduce `metacheck_matrix_x`, `metacheck_matrix_z`, and `metacheck_matrix`.
- Introduce `generator_polynomial`.
- Move the following codes from `QuantumClifford.ECC` to `QECCore`: `ReedMuller`, `RecursiveReedMuller`, `QuantumReedMuller`, `Hamming`, `Golay`, `Triangular488 `, `Triangular666 `, `Gottesman`.
- Add classical `Goppa`
- Add classical Gallager's LDPC code
- Add cyclic quantum Tanner graph product codes
- Add `[[n² + m²,(n - rank([C ∣ M]))² + (m − rank([C ∣ M]ᵀ))², d]]` quantum Tillich-Zémor `random_TillichZemor_code`
- Add `Delfosse-Reichardt` codes from classical self-orthogonal `Reed-Muller` seed codes
- Add `[[4p, 2(p − 2), 4]]` Delfosse-Reichardt repetition `DelfosseReichardtRep`
- Add `[[8p, 4p − 2, 3]]` Delfosse-Reichardt Generalized `[[8,2,3]]` `DelfosseReichardt823`
- In an Oscar extension for `QECCore`
    - Add `BivariateBicycleViaCirculantMat`
- Create a Nemo.jl package extension for binary matrix `rank` computation.

## v0.1.1

- Add `AbstractDistanceAlg` for dispatch for `distance`

## v0.1.0

- First release, moving basic definitions from QuantumClifford.jl
