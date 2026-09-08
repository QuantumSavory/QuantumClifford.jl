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
    The concrete parity matrices label ``x^i y^j`` by the canonical index
    ``i + \\ell j + 1``, for ``0 \\leq i < \\ell`` and ``0 \\leq j < m``.
    Thus, equal constructor arguments always give the same check and qubit labels.

Here is an example of the `[[80, 2, 8]]` ZSZ code from Table I of [guo2025zsz](@cite).
The parameters `l=5`, `m=8`, and `q=2` correspond to the semidirect product presentation ``\\langle x, y \\mid x^5=1, y^8=1, y x y^{-1} = x^2 \\rangle``.
Note that ZSZ codes are generally asymmetric, meaning the X-distance and Z-distance differ (`d_X != d_Z`). For instance, the `[[80, 2, 8]]` code below has `d_X = 10` and `d_Z = 8`.

```jldoctest
julia> using Oscar, QuantumClifford.ECC;

julia> import HiGHS;

julia> c = ZSZ(5, 8, 2, [(0,0),(4,4),(4,1)], [(0,0),(3,0),(2,7)]);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(80, 2, 8)
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
        powermod(q, m, l) == 1 || throw(ArgumentError("Condition q^m ≡ 1 (mod l) not satisfied for l=$l, m=$m, q=$q"))
        new(l, m, q, A, B)
    end
end

function parity_matrix_xz(c::ZSZ)
    l, m = c.l, c.m
    group_order = l * m
    Hx = zeros(Bool, group_order, 2 * group_order)
    Hz = zeros(Bool, group_order, 2 * group_order)
    A = [(mod(i, l), mod(j, m)) for (i, j) in c.A]
    B = [(mod(i, l), mod(j, m)) for (i, j) in c.B]
    q_powers = [powermod(c.q, j, l) for j in 0:m-1]

    # Elements use the normal form x^i*y^j and index 1 + i + l*j.
    # Their product is (i, j)*(k, r) = (i + q^j*k mod l, j + r mod m).
    # The LPCode convention is Hx = [R(A) L(B)] and Hz = [L(B)' R(A)'].
    for j in 0:m-1
        qj = q_powers[j + 1]
        for i in 0:l-1
            row = 1 + i + l * j
            for (ai, aj) in A
                col = 1 + mod(i + qj * ai, l) + l * mod(j + aj, m)
                Hx[col, row] ⊻= true
                Hz[row, group_order + col] ⊻= true
            end
            for (bi, bj) in B
                col = 1 + mod(bi + q_powers[bj + 1] * i, l) + l * mod(bj + j, m)
                Hx[col, group_order + row] ⊻= true
                Hz[row, col] ⊻= true
            end
        end
    end
    return Hx, Hz
end

parity_matrix_x(c::ZSZ) = first(parity_matrix_xz(c))
parity_matrix_z(c::ZSZ) = last(parity_matrix_xz(c))

function parity_matrix(c::ZSZ)
    return parity_matrix(CSS(parity_matrix_xz(c)...))
end

code_n(c::ZSZ) = 2 * c.l * c.m
