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

!!! note
    This function is simply a convenience wrapper that handles argument conversions before calling [`two_block_group_algebra_code`](@ref). Notably, it uses Cayley's theorem to compute a group isomorphism of a finitely presented group as a permutation group, enabling the efficient construction of larger blocklength ZSZ codes from Table I [guo2025zsz](@cite).

Here is an example of the `[[80, 2]]` ZSZ code from Table I of [guo2025zsz](@cite):

```jldoctest
julia> using Oscar, QuantumClifford.ECC;

julia> c = ZSZ(5, 8, 2, [(0,0),(4,4),(4,1)], [(0,0),(3,0),(2,7)]);

julia> code_n(c), code_k(c)
(80, 2)
```

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

function zsz_to_lpcode(c::ZSZ)
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
    return two_block_group_algebra_code(a, b)
end

parity_matrix_x(c::ZSZ) = parity_matrix_x(zsz_to_lpcode(c))
parity_matrix_z(c::ZSZ) = parity_matrix_z(zsz_to_lpcode(c))

# override because the generic method calls parity_matrix_x and parity_matrix_z
# separately, each building a new lpcode with a different non-deterministic
# isomorphism from Oscar.isomorphism(PermGroup, Q). building once ensures
# Hx and Hz share the same group element ordering so CSS commutativity holds.
function parity_matrix(c::ZSZ)
    lp = zsz_to_lpcode(c)
    return parity_matrix(CSS(parity_matrix_x(lp), parity_matrix_z(lp)))
end

code_n(c::ZSZ) = 2 * c.l * c.m
