"""
# Direct Product of Groups

The direct product of groups is instrumental in constructing group algebra of two-block
group algebra code. Lin and Pryadko illustrate this method in Appendix C, Table 2 of
[lin2024quantum](@cite), where they utilize the direct product of two cyclic groups,
expressed as `C₂ₘ = Cₘ × C₂`, with an order of `2m`.

`Hecke.jl` contains only abelian groups and a list of all finite groups of order up to 100.
`Oscar.jl` brings in comprehensive functionality for computational group theory, including
support for **arbitrary finitely presented groups** (groups of the form `⟨X | S⟩`. `Oscar.jl`
supports the **direct product** operation between two or more arbitrary **general** groups,
including non-abelian groups such as `alternating_group`, `dihedral_group`, `symmetric_group`,
and even arbitrary finitely presented groups (e.g., `free_group`). This capability is not
available in `Hecke.jl`. The 2BGA codes discovered in [lin2024quantum](@cite) rely on direct
products of two or more *general* groups, which necessitate the use of `Oscar.direct_product`.

The schematic below illustrates the limitations of `Hecke.direct_product` compared to
`Oscar.direct_product`:

```@raw html
<div class="mermaid">
graph TB
    root[Direct Product of Groups]

    root --> A[Hecke.direct_product]
    root --> B[Oscar.direct_product]

    %% Hecke Branch
    A --> A1[Supports mostly abelian groups and list of finite groups]
    A1--> A2[abelian_group symmetric_group small_group]
    A2 --> A3[C × C, C × S]

    %% Oscar Branch
    B --> B1[Supports finite general groups, including non-abelian groups]
    B1--> B2[alternating_group <br> dihedral_group <br> free_group <br> cyclic_group <br> permutation_group <br>  quaternion_group <br> symmetric_group <br> abelian_group<br> small_group]
    B2--> B3[A × C,  A × D <br> D × C, D × D <br> F × F, F × A, F × D <br> F × S, F × C <br> C × C, C × S]
</div>
```

# Example

The [[56, 28, 2]] abelian 2BGA code from Appendix C, Table II in [lin2024quantum](@cite)
can be constructed using the direct product of two cyclic groups. Specifically, the group
`C₂₈` of order `l = 28` can be represented as `C₁₄ × C₂`, where the first group has order
`m = 14` and the second group has order `n = 2`.

```jldoctest directprod
julia> import Oscar: cyclic_group, small_group_identification, describe, order; # hide

julia> import Hecke: gens, quo, group_algebra, GF, one, direct_product, sub; # hide

julia> m = 14; n = 2;

julia> C₁₄ = cyclic_group(m);

julia> C₂ = cyclic_group(n);

julia> G = direct_product(C₁₄, C₂);

julia> GA = group_algebra(GF(2), G);

julia> x, s = gens(GA)[1], gens(GA)[3];

julia> a = [one(GA), x^7];

julia> b = [one(GA), x^7, s, x^8, s * x^7, x];

julia> c = twobga_from_direct_product(a, b, GA);

julia> order(G)
28

julia> code_n(c), code_k(c)
(56, 28)

julia> describe(G), small_group_identification(G)
("C14 x C2", (28, 4))
```

!!! note When using the direct product of two cyclic groups, it is essential to verify
the group presentation `Cₘ = ⟨x, s | xᵐ = s² = xsx⁻¹s⁻¹ = 1⟩` is satisfied, where the
order is `2m`. Ensure that the selected generators have the correct orders of `m = 14`
and `n = 2`, respectively. If the group presentation is not satisfied, the resulting
group algebra over `GF(2)` will not represent the intended group, `C₂₈ = C₁₄ × C₂`. In
addition, `Oscar.sub` can be used to determine if `H` is a subgroup of `G` and to
confirm that both `C₁₄` and `C₂` are subgroups of `C₂₈`.

```jldoctest directprod
julia> order(gens(G)[1])
14

julia> order(gens(G)[3])
2

julia> x^14 == s^2 == x * s * x^-1 * s^-1
true

julia> H, _  = sub(G, [gens(G)[1], gens(G)[3]]);

julia> H == G
true
```
"""
function twobga_from_direct_product(a_elts::VectorDirectProductGroupElem, b_elts::VectorDirectProductGroupElem, F2G::DirectProductGroupAlgebra)
    a = sum(F2G(x) for x in a_elts)
    b = sum(F2G(x) for x in b_elts)
    c = two_block_group_algebra_codes(a,b)
    return c
end
