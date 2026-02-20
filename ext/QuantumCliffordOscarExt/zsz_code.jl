"""
    $TYPEDEF

A zero-shot ZSZ code from [guo2025zsz](@cite).
Returns a CSS code constructed from the two-block group algebra over the semidirect product group ``\\mathbb{Z}_\\ell \\rtimes_q \\mathbb{Z}_m``.

### Fields
    $TYPEDFIELDS
"""
struct ZSZCode <: AbstractCSSCode
    """Order of the normal cyclic group ``\\mathbb{Z}_\\ell``"""
    l::Int
    """Order of the acting cyclic group ``\\mathbb{Z}_m``"""
    m::Int
    """The action parameter `q` such that ``y x y^{-1} = x^q``"""
    q::Int
    """First polynomial A represented as a list of `(i, j)` exponent tuples for ``x^i y^j``"""
    A::Vector{Tuple{Int, Int}}
    """Second polynomial B represented as a list of `(i, j)` exponent tuples for ``x^i y^j``"""
    B::Vector{Tuple{Int, Int}}

    function ZSZCode(l::Int, m::Int, q::Int, A::Vector{Tuple{Int,Int}}, B::Vector{Tuple{Int,Int}})
        if powermod(q, m, l) != 1
            throw(ArgumentError("Condition q^m = 1 (mod l) not satisfied for l=$l, m=$m, q=$q"))
        end
        new(l, m, q, A, B)
    end
end

function _zsz_to_css(c::ZSZCode)
    G = free_group(2)
    x, y = gens(G)
    rels = [x^c.l, y^c.m, y*x*y^-1 * x^-c.q]
    Q, _ = quo(G, rels)
    F2G = group_algebra(GF(2), Q)
    qx, qy = gens(Q)

    A_elts = typeof(qx)[]
    for (i, j) in c.A
        push!(A_elts, qx^i * qy^j)
    end
    B_elts = typeof(qx)[]
    for (i, j) in c.B
        push!(B_elts, qx^i * qy^j)
    end
    
    return twobga_from_fp_group(A_elts, B_elts, F2G)
end

parity_matrix_x(c::ZSZCode) = parity_matrix_x(_zsz_to_css(c))

parity_matrix_z(c::ZSZCode) = parity_matrix_z(_zsz_to_css(c))

code_n(c::ZSZCode) = 2 * c.l * c.m
