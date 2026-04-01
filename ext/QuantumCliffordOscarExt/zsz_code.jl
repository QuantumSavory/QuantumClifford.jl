"""
    $TYPEDEF

ZSZ codes are single-shot decodable [`two_block_group_algebra_code`](@ref)s built from the semidirect product of groups ``\\mathbb{Z}_\\ell \\rtimes_q \\mathbb{Z}_m`` [guo2025zsz](@cite).

This code is defined by the group presentation:

```math
\\begin{aligned}
\\langle x, y \\mid x^\\ell = 1, y^m = 1, y x y^{-1} = x^q \\rangle
\\end{aligned}
```

Notably, it is an instance of a [`two_block_group_algebra_code`](@ref) code with this specific presentation. While it lacks explicit *metachecks*, it exhibits single-shot properties (e.g., self-correction with passive greedy decoding) due to strong error confinement stemming from small-set expansion in its Tanner graph [guo2025zsz](@cite).

This particular function is nothing more than a simple wrapper that takes care of argument conversions for [`two_block_group_algebra_code`](@ref).
Of note, the polynomials here are given as lists of `(i, j)` exponent tuples for ``x^i y^j``.

### Fields
    $TYPEDFIELDS
"""
struct ZSZ <: AbstractCSSCode
    """Order of the cyclic group ``\\mathbb{Z}_\\ell``"""
    l::Int
    """Order of the cyclic group ``\\mathbb{Z}_m``"""
    m::Int
    """The parameter `q` such that ``y x y^{-1} = x^q``"""
    q::Int
    """First polynomial A represented as a list of `(i, j)` exponent tuples for ``x^i y^j``"""
    A::Vector{Tuple{Int, Int}}
    """Second polynomial B represented as a list of `(i, j)` exponent tuples for ``x^i y^j``"""
    B::Vector{Tuple{Int, Int}}

    function ZSZ(l::Int, m::Int, q::Int, A::Vector{Tuple{Int,Int}}, B::Vector{Tuple{Int,Int}})
        if powermod(q, m, l) != 1
            throw(ArgumentError("Condition q^m = 1 (mod l) not satisfied for l=$l, m=$m, q=$q"))
        end
        new(l, m, q, A, B)
    end
end

"""Constructs parity check matrices for the ZSZ code.

Converts the finitely presented group to a permutation group before constructing
the group algebra, because `group_algebra(GF(2), FPGroup)` hangs for larger group
orders while `PermGroup` works efficiently."""
function parity_matrix_xz(c::ZSZ)
    G = free_group(2)
    x, y = gens(G)
    rels = [x^c.l, y^c.m, y*x*y^-1 * x^-c.q]
    Q, _ = quo(G, rels)
    iso = isomorphism(PermGroup, Q)
    Qp = codomain(iso)
    F2G = group_algebra(GF(2), Qp)
    qx, qy = gens(Q)
    a = sum(F2G(iso(qx^i * qy^j)) for (i, j) in c.A)
    b = sum(F2G(iso(qx^i * qy^j)) for (i, j) in c.B)
    c2 = two_block_group_algebra_code(a, b)
    return parity_matrix_xz(c2)
end

parity_matrix_x(c::ZSZ) = parity_matrix_xz(c)[1]
parity_matrix_z(c::ZSZ) = parity_matrix_xz(c)[2]

"""Need to override `parity_matrix` because the generic method calls `parity_matrix_x`
and `parity_matrix_z` separately, each of which would build a new group algebra with
a different isomorphism. This way it is built once."""
function parity_matrix(c::ZSZ)
    hx, hz = parity_matrix_xz(c)
    return parity_matrix(CSS(hx, hz))
end

code_n(c::ZSZ) = 2 * c.l * c.m
