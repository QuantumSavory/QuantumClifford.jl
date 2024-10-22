"""
# Specific Group Presentations

For quantum error-correcting codes, such as the two-block algebra (2BGA) code, functionalities
like designing specific group presentations for both abelian and non-abelian groups are crucial.
These specific presentations are the key ingredient for construction of group algebra of 2BGA
with abelian and non-abelian groups.

# Example

[[96, 12, 10]] 2BGA code from Table I of [lin2024quantum](@cite) with group presentation
`⟨r, s|s⁶ = r⁸ = r⁻¹srs = 1⟩` and group structure C₂ × (C₃ ⋉ C₈).

```jldoctest finitegrp
julia> import Oscar: free_group, small_group_identification, describe; # hide

julia> import Hecke: gens, quo, group_algebra, GF, one; # hide

julia> F = free_group(["r", "s"]);

julia> r, s = gens(F); # generators

julia> G, = quo(F, [s^6, r^8, r^(-1) * s * r * s]);  # relations

julia> GA = group_algebra(GF(2), G);

julia> r, s = gens(G);

julia> a = [one(G), r, s^3 * r^2, s^2 * r^3];

julia> b = [one(G), r, s^4 * r^6, s^5 * r^3];

julia> c = twobga_from_fp_group(a, b, GA);

julia> code_n(c), code_k(c)
(96, 12)

julia> describe(G), small_group_identification(G)
("C2 x (C3 : C8)", (48, 9))
```

# Cyclic Groups

Cyclic groups with specific group presentations `Cₘ = ⟨x, s|xᵐ = s² = xsx⁻¹s⁻¹ = 1⟩` where order
is `2m` are supported.

To construct a group algebra for abelian cyclic groups, specify the group's presentation `⟨S|R⟩`
using its generators `S` and defining relations `R`.

# Example

[[56, 28, 2]] 2BGA code from Appendix C, Table II of [lin2024quantum](@cite) for abelian
cyclic group `C₂₈`.

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

julia> code_n(c), code_k(c)
(56, 28)
```

# Dihedral Groups

Dihedral groups with specific group presentations `Dₘ = ⟨r, s|rᵐ = s² = (rs)² = 1⟩` where order
is `2m` are supported.

To construct a group algebra for non-abelian Dihedral groups, specify the group's presentation
`⟨S|R⟩` using its generators `S` and defining relations `R`.

# Example

[[24, 8, 3]] 2BGA code from Appendix C, Table III of [lin2024quantum](@cite) for non-abelian
Dihedral group `D₆`.

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

julia> code_n(c), code_k(c)
(24, 8)
```

"""
function twobga_from_fp_group(a_elts::Vector{FPGroupElem}, b_elts::Vector{FPGroupElem}, F2G::GroupAlgebra{FqFieldElem, FPGroup, FPGroupElem})
    a = sum(F2G(x) for x in a_elts)
    b = sum(F2G(x) for x in b_elts)
    c = two_block_group_algebra_codes(a,b)
    return c
end
