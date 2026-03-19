"""
    $TYPEDEF

A single-shot ZSZ code
Returns a CSS code constructed from the two-block group algebra (2BGA) over the semidirect product group ``\\mathbb{Z}_\\ell \\rtimes_q \\mathbb{Z}_m``.

This code is defined by the group presentation:
``\\langle x, y \\mid x^\\ell = 1, y^m = 1, y x y^{-1} = x^q \\rangle``

It is an instance of a 2BGA code with this specific presentation. While it lacks an extensive number of redundant parity checks (metachecks), the paper claims it exhibits single-shot properties (e.g., self-correction with passive greedy decoding) due to strong error confinement stemming from small-set expansion in its Tanner graph.

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

# builds an LPCode from the zsz parameters.
# convert the finitely presented group to a permutation group before
# constructing the group algebra, because group_algebra(GF(2), FPGroup)
# hangs for larger group orders while PermGroup works fine.
function _zsz_to_lpcode(c::ZSZ)
    G = free_group(2)
    x, y = gens(G)
    rels = [x^c.l, y^c.m, y*x*y^-1 * x^-c.q]
    Q, _ = quo(G, rels)
    # fp group -> perm group so that group_algebra doesnt hang
    iso = Oscar.isomorphism(PermGroup, Q)
    Qp = Oscar.codomain(iso)
    F2G = group_algebra(GF(2), Qp)
    qx, qy = gens(Q)
    a = sum(F2G(iso(qx^i * qy^j)) for (i, j) in c.A)
    b = sum(F2G(iso(qx^i * qy^j)) for (i, j) in c.B)
    return two_block_group_algebra_code(a, b)
end

parity_matrix_xz(c::ZSZ) = parity_matrix_xz(_zsz_to_lpcode(c))
parity_matrix_x(c::ZSZ) = parity_matrix_xz(c)[1]
parity_matrix_z(c::ZSZ) = parity_matrix_xz(c)[2]
# need to override parity_matrix here because the generic method
# calls parity_matrix_x and parity_matrix_z separately, and each of
# those builds a new lpcode with a different isomorphism.
# this way it is built once and both hx, hz are obtained from the same construction.
function parity_matrix(c::ZSZ)
    hx, hz = parity_matrix_xz(c)
    return parity_matrix(CSS(hx, hz))
end

code_n(c::ZSZ) = 2 * c.l * c.m
