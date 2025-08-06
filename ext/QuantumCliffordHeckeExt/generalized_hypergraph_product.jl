"""
    $TYPEDEF

A Generalized Hypergraph Product (GHP) CSS code introduced in [panteleev2021degenerate](@cite).

The GHP code is constructed from a matrix ``A \\in M_{m \\times n}(R)``, where
``R \\subseteq M_\\ell(\\mathbb{F}_2)`` is a ring of binary circulant matrices,
and a binary matrix ``b \\in M_\\ell(\\mathbb{F}_2)`` such that every element of
``R`` commutes with ``b``. The code is defined by the block parity-check matrices:

```math
\\begin{aligned}
    H_X = \\begin{bmatrix} A & b \\cdot I_m \\end{bmatrix}, \\quad
    H_Z = \\begin{bmatrix} b^T \\cdot I_n & A^* \\end{bmatrix}
\\end{aligned}
```

where ``A^*`` is the transpose of ``A`` with each entry polynomial-reversed modulo
``x^\\ell - 1``, and ``I_m, I_n`` are identity matrices over ``R``.

Here is an example of `[[882, 24, 18 ≤ d ≤ 24]]` code from Appendix B of [panteleev2021degenerate](@cite).

```jldoctest
julia> import Hecke: polynomial_ring, GF, quo, matrix; using QuantumClifford.ECC;

julia> R, x = polynomial_ring(GF(2), :x);

julia> l = 63; n = 7;

julia> S, _ =  quo(R, x^l - 1);

julia> A = matrix(S, n, n,
           [x^27  0     0     0     0     1     x^54
            x^54  x^27  0     0     0     0     1
            1     x^54  x^27  0     0     0     0
            0     1     x^54  x^27  0     0     0
            0     0     1     x^54  x^27  0     0
            0     0     0     1     x^54  x^27  0
            0     0     0     0     1     x^54  x^27]);

julia> b = S(1 + x + x^6);

julia> c = GeneralizedHyperGraphProductCode(A, b, l);

julia> code_n(c), code_k(c)
(882, 24)
```

Here is an example of `[[1270, 28, 16 ≤ d ≤ 46]]` code from Appendix B of [panteleev2021degenerate](@cite).

```jldoctest
julia> import Hecke: polynomial_ring, GF, quo, matrix; using QuantumClifford.ECC;

julia> R, x = polynomial_ring(GF(2), :x);

julia> l = 127; n = 5;

julia> S, _ =  quo(R, x^l - 1);

julia> A = matrix(S, n, n,
           [1     0     x^51  x^52  0
            0     1     0     x^111 x^20
            1     0     x^98  0     x^122
            1     x^80  0     x^119 0 
            0     1     x^5   0     x^106]);

julia> b = S(1 + x + x^7);

julia> c = GeneralizedHyperGraphProductCode(A, b, l);

julia> code_n(c), code_k(c)
(1270, 28)
```

### Fields
    $TYPEDFIELDS
"""
struct GeneralizedHyperGraphProductCode <: AbstractCSSCode
    """The matrix ``A \\in M_{n \\times n}(R)``, where ``R = \\mathbb{F}_2[x]/(x^\\ell - 1)``. Each entry in `A` represents a polynomial modulo ``x^\\ell - 1``, defining a circulant block."""
    A::MatSpaceElem{EuclideanRingResidueRingElem{FqPolyRingElem}}
    """The polynomial ``b(x) \\in R``, generating a binary circulant matrix that commutes with all elements of `R`. """
    b::EuclideanRingResidueRingElem{FqPolyRingElem}
    """The number of rows and columns in each binary circulant block."""
    l::Int

    function GeneralizedHyperGraphProductCode(A, b, l)
        size(A, 1) == size(A, 2) || throw(ArgumentError("A must be square"))
        parent(A[1,1]) == parent(b) || throw(ArgumentError("A and b must be over the same ring"))
        new(A, b, l)
    end
end

function _poly_transpose(p, l)
    plift = lift(p)
    deg_p = degree(plift)
    rev = zero(base_ring(parent(p)))
    x = gen(parent(p))
    for i in 0:deg_p
        rev += coeff(plift, i)*x^(mod(-i, l))
    end
    return parent(p)(rev)
end

_bᵢ(b, n) = matrix(parent(b), n, n, [i == j ? b : zero(parent(b)) for i in 1:n, j in 1:n])

function _circulant_matrix(poly, l)
    coeffs = [coeff(poly, i) for i in 0:l-1]
    return matrix(GF(2), l, l, [coeffs[mod1(j-i+1, l)] for i in 1:l, j in 1:l])
end

function _polynomial_matrix_to_circulant_matrix(H_poly, l)
    m, n = size(H_poly)
    H = zero_matrix(GF(2), m*l, n*l)
    for i in 1:m, j in 1:n
        poly = lift(H_poly[i, j])
        H_block = _circulant_matrix(poly, l)
        H[(i-1)*l+1:i*l, (j-1)*l+1:j*l] = H_block
    end
    H = [Int(lift(ZZ, H[i,j])) for i in 1:nrows(H), j in 1:ncols(H)]
    return H
end

function parity_matrix_xz(c::GeneralizedHyperGraphProductCode)
    n = size(c.A, 1)
    Aᵗʳ = matrix(parent(c.b), n, n, [_poly_transpose(c.A[j, i], c.l) for i in 1:n, j in 1:n])
    bᵢ = _bᵢ(c.b, n)
    bᵢᵀ = _bᵢ(_poly_transpose(c.b, c.l), n)
    hx = hcat(c.A, bᵢ)
    hz = hcat(bᵢᵀ, Aᵗʳ)
    hx = _polynomial_matrix_to_circulant_matrix(hx, c.l)
    hz = _polynomial_matrix_to_circulant_matrix(hz, c.l)
    return hx, hz
end

parity_matrix_x(c::GeneralizedHyperGraphProductCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::GeneralizedHyperGraphProductCode) = parity_matrix_xz(c)[2]
