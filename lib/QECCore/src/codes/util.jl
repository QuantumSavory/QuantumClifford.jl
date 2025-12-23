"""
    $TYPEDEF

Represents the CSS quantum code `Q(G₁ × G₂)` constructed from two binary codes with
parity-check matrices `H₁` and `H₂`, using the hypergraph product formulation introduced
by [tillich2013quantum](@cite).

This construction corresponds to a specific product of Tanner graphs:
- Let `G₁ = T(V₁, C₁, E₁)` and `G₂ = T(V₂, C₂, E₂)` be the Tanner graphs of `H₁` and `H₂`.
- The product graph `G₁ × G₂` has vertex set `(V₁ × V₂) ∪ (C₁ × C₂)` and check set `(C₁ × V₂) ∪ (V₁ × C₂)`.
- The Tanner subgraphs `G₁ ×ₓ G₂` and `G₁ ×𝓏 G₂` define classical codes `Cₓ` and `C𝓏` used in the CSS construction.

The `hgp(H₁, H₂)` function algebraically realizes this graph-theoretic product using Kronecker operations,
yielding the `X`- and `Z`-type parity-check matrices:

- `H_X = [H₁ ⊗ I  |  I ⊗ H₂ᵗ]` corresponds to `G₁ ×ₓ G₂`
- `H_Z = [I ⊗ H₂  |  H₁ᵗ ⊗ I]` corresponds to `G₁ ×𝓏 G₂`

These matrices ensure `H_X * H_Zᵗ = 0`, satisfying the CSS condition.

See: [tillich2013quantum](@cite), Section 4.3 — “The hypergraph connection, product codes”

# 𝑄(𝐺₁ × 𝐺₂)

The `𝑄(𝐺₁ × 𝐺₂)` quantum LDPC codes represent a broader generalization of **quantum expander codes**
which are derived from the Leverrier-Tillich-Zémor construction [tillich2013quantum](@cite).

```jldoctest examples
julia> using QuantumClifford; using QuantumClifford.ECC; using QECCore

julia> H1 = [1 0 1 0; 0 1 0 1; 1 1 0 0];

julia> H2 = [1 1 0; 0 1 1];

julia> c = parity_checks(hgp(H1, H2)...)
+ X_____X_____X_____
+ _X_____X____XX____
+ __X_____X____X____
+ ___X_____X____X___
+ ____X_____X___XX__
+ _____X_____X___X__
+ X__X____________X_
+ _X__X___________XX
+ __X__X___________X
+ ZZ__________Z___Z_
+ _ZZ__________Z___Z
+ ___ZZ_________Z_Z_
+ ____ZZ_________Z_Z
+ ______ZZ____Z_____
+ _______ZZ____Z____
+ _________ZZ___Z___
+ __________ZZ___Z__

julia>  code_n(c), code_k(c)
(18, 1)
```

# Quantum Expander code

The `𝑄(𝐺₁ × 𝐺₂)` code is more general than the **standard quantum expander code**
[leverrier2015quantum](@cite) construction. The quantum expander code construction
corresponds to the specific case where `G = G1 = G2`​.

```jldoctest examples
julia> H = parity_matrix(RepCode(3));

julia> c = parity_checks(hgp(H, H)...)
+ X__X_____X_X______
+ _X__X____XX_______
+ __X__X____XX______
+ ___X__X_____X_X___
+ ____X__X____XX____
+ _____X__X____XX___
+ X_____X________X_X
+ _X_____X_______XX_
+ __X_____X_______XX
+ ZZ_______Z_____Z__
+ _ZZ_______Z_____Z_
+ Z_Z________Z_____Z
+ ___ZZ____Z__Z_____
+ ____ZZ____Z__Z____
+ ___Z_Z_____Z__Z___
+ ______ZZ____Z__Z__
+ _______ZZ____Z__Z_
+ ______Z_Z_____Z__Z

julia>  code_n(c), code_k(c)
(18, 2)
```

See also: [`CyclicQuantumTannerGraphProduct`](@ref)
"""
function hgp(h₁,h₂)
    r₁, n₁ = size(h₁)
    r₂, n₂ = size(h₂)
    hx = hcat(kron(h₁, LinearAlgebra.I(n₂)), kron(LinearAlgebra.I(r₁), h₂'))
    hz = hcat(kron(LinearAlgebra.I(n₁), h₂), kron(h₁', LinearAlgebra.I(r₂)))
    hx, hz
end
