"""
$TYPEDSIGNATURES

Constructing two block group algebra codes by specifying the direct product to be used.
See also the more general [`two_block_group_algebra_codes`](@ref).

Two block group algebra codes are constructed by choosing a group (and specific generators),
then choosing two polynomials made out of these generators,
then piping these two polynomials as the elements of `1×1` matrices to the
lifted product code constructors.

The Hecke library, for which we already have an extension, provides for a fairly easy way
to construct such polynomials for many abelian and small groups.
See [`two_block_group_algebra_codes`](@ref) for those capabilities.

However, more esoteric groups can be specified as the direct product of other groups.
To support arbitrary direct products we use Oscar, which builds upon Hecke.
Oscar supports the **direct product** operation between two or more arbitrary **general** groups,
including non-abelian groups such as `alternating_group`, `dihedral_group`, `symmetric_group`,
and even arbitrary finitely presented groups (e.g., `free_group`). This capability is not
available in `Hecke.jl`. The 2BGA codes discovered in [lin2024quantum](@cite) rely on direct
products of two or more *general* groups, which necessitate the use of `Oscar.direct_product`.

This particular function is nothing more than a simple wrapper that takes care of argument conversions.
Of note, the polynomials here are given as lists of monomials.

Of course, if you are comfortable with Oscar, you can use [`two_block_group_algebra_codes`](@ref) directly.

See also: [`two_block_group_algebra_codes`](@ref), [`twobga_from_fp_group`](@ref)

## Examples

The [[56, 28, 2]] abelian 2BGA code from Appendix C, Table II in [lin2024quantum](@cite)
can be constructed using the direct product of two cyclic groups. Specifically, the group
`C₂₈` of order `l = 28` can be represented as `C₁₄ × C₂`, where the first group has order
`m = 14` and the second group has order `n = 2`.

```jldoctest directprod
julia> import Oscar: cyclic_group, small_group_identification, describe, order

julia> import Hecke: gens, quo, group_algebra, GF, one, direct_product, sub

julia> using QuantumClifford, QuantumClifford.ECC

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

!!! danger
    When using the direct product, there isn't necessarily a unique set of generators.
    It is essential to verify that Oscar is providing you with the generators you expect,
    e.g. for a cycling group that you have the presentation `Cₘ = ⟨x, s | xᵐ = s² = xsx⁻¹s⁻¹ = 1⟩`.
    For situations where the generators provided by Oscar are not the ones you want,
    you can also use [`twobga_from_fp_group`](@ref) where you specify the group presentation directly.

As a verification that you have the correct generators,
`Oscar.sub` can be used to determine if `H` is a subgroup of `G` and to
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
