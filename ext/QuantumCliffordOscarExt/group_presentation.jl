"""
# Specific Group Presentations

For quantum error-correcting codes like the two-block group algebra (2BGA) code, designing specific
group presentations for both abelian and non-abelian groups is crucial. These presentations are essential
for constructing the group algebra of the 2BGA code for a given finite general group, `G`.

Lin and Pryadko, in their seminal paper titled "Quantum Two-Block Group Algebra Codes" [lin2024quantum](@cite),
employ specific presentations that necessitate the use of `Oscar.free_group`. The diagram below distinguishes
between small groups (`Hecke/Oscar.small_group`) and finitely presented groups (`Oscar.free_group`) by highlighting
the existence of extra relations in their presentations.

```@raw html
<div class="mermaid">
graph TD
    A[Group Presentation ⟨S ∣ R⟩] --> B{Are there <br> extra relations?}
    B -- No --> C[Small groups <br> Hecke/Oscar.small_group]
    C --> D[Independent generators]
    C --> E["Example: <br> ⟨r, s ∣ s⁴, r⁹⟩"]
    B -- Yes --> F[Finitely presented groups <br> Oscar.free_group]
    F --> G[Defined by interactions]
    F --> H["Example: <br> ⟨r, s ∣ s⁴, r⁹, s⁻¹rsr⟩"]
</div>
```

# Example

The [[96, 12, 10]] 2BGA code from Table I in [lin2024quantum](@cite) has the group presentation
`⟨r, s | s⁶ = r⁸ = r⁻¹srs = 1⟩` and a group structure of `C₂ × (C₃ ⋉ C₈)`.

```jldoctest finitegrp
julia> import Oscar: free_group, small_group_identification, describe, order; # hide

julia> import Hecke: gens, quo, group_algebra, GF, one; # hide

julia> F = free_group(["r", "s"]);

julia> r, s = gens(F); # generators

julia> G, = quo(F, [s^6, r^8, r^(-1) * s * r * s]);  # relations

julia> GA = group_algebra(GF(2), G);

julia> r, s = gens(G);

julia> a = [one(G), r, s^3 * r^2, s^2 * r^3];

julia> b = [one(G), r, s^4 * r^6, s^5 * r^3];

julia> c = twobga_from_fp_group(a, b, GA);

julia> order(G)
48

julia> code_n(c), code_k(c)
(96, 12)

julia> describe(G), small_group_identification(G)
("C2 x (C3 : C8)", (48, 9))
```

# Cyclic Groups

Cyclic groups with specific group presentations, given by `Cₘ = ⟨x, s | xᵐ = s² = xsx⁻¹s⁻¹ = 1⟩`,
where the order is `2m`, are supported.

To construct the group algebra for a cyclic group, specify the group presentation `⟨S | R⟩`, using
its generators `S` and defining relations `R`.

# Example

The [[56, 28, 2]] abelian 2BGA code from Appendix C, Table II in [lin2024quantum](@cite) is constructed using
the cyclic group `C₂₈ = C₁₄ × C₂`.

```jldoctest finitegrp
julia> m = 14;

julia> F = free_group(["x", "s"]);

julia> x, s = gens(F); # generators

julia> G, = quo(F, [x^m, s^2, x * s * x^-1 * s^-1]); # relations

julia> GA = group_algebra(GF(2), G);

julia> x, s = gens(G);

julia> a = [one(G), x^7];

julia> b = [one(G), x^7, s, x^8, s * x^7, x];

julia> c = twobga_from_fp_group(a, b, GA);

julia> order(G)
28

julia> code_n(c), code_k(c)
(56, 28)

julia> describe(G), small_group_identification(G)
("C14 x C2", (28, 4))
```

# Dihedral Groups

Dihedral groups with specific group presentations, given by `Dₘ = ⟨r, s | rᵐ = s² = (rs)² = 1⟩`,
where the order is `2m`, are supported.

To construct the group algebra for a dihedral group, specify the group presentation `⟨S | R⟩`
using its generators `S` and defining relations `R`.

# Example

The [[24, 8, 3]] 2BGA code from Appendix C, Table III in [lin2024quantum](@cite) is constructed
using the dihedral group `D₆ = C₆ ⋉ C₂`.

```jldoctest finitegrp
julia> m = 6;

julia> F = free_group(["r", "s"]);

julia> r, s = gens(F); # generators

julia> G, = quo(F, [r^m, s^2, (r*s)^2]); # relations

julia> GA = group_algebra(GF(2), G);

julia> r, s = gens(G);

julia> a = [one(G), r^4];

julia> b = [one(G), s*r^4, r^3, r^4, s*r^2, r];

julia> c = twobga_from_fp_group(a, b, GA);

julia> order(G)
12

julia> code_n(c), code_k(c)
(24, 8)

julia> describe(G), small_group_identification(G)
("D12", (12, 4))
```
"""
function twobga_from_fp_group(a_elts::VectorFPGroupElem, b_elts::VectorFPGroupElem, F2G::FqFieldFPGroupAlgebra)
    a = sum(F2G(x) for x in a_elts)
    b = sum(F2G(x) for x in b_elts)
    c = two_block_group_algebra_codes(a,b)
    return c
end
