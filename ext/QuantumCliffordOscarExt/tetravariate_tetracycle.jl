"""
    $TYPEDEF

A novel generalization of [`TrivariateTricycleCode`](@ref) built from
tetravariate polynomials A, B, C, D in the multivariate polynomial quotient
ring ``\\frac{\\mathbb{F}_2[w, x, y, z]}{\\langle w^\\ell-1, x^m-1, y^p-1 \\rangle, y^p-1 \\rangle}``.

Here is the novel `[[96, 6, d]]` tetravariate tetracycle code.

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p, r = 2, 2, 2, 2;

julia> R, (w, x, y, z) = polynomial_ring(GF(2), [:w, :x, :y, :z]);

julia> I = ideal(R, [w^l - 1, x^m - 1, y^p - 1, z^r - 1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + w);

julia> B = S(1 + x);

julia> C = S(1 + y);

julia> D = S(1 + z);

julia> c = TetravariateTetracycleCode(l, m, p, r, A, B, C, D);

julia> code_n(c), code_k(c)
(96, 6)
```

Here is the novel `[[486, 12, d]]` tetravariate tetracycle code.

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p, r = 3, 3, 3, 3;

julia> R, (w, x, y, z) = polynomial_ring(GF(2), [:w, :x, :y, :z]);

julia> I = ideal(R, [w^l - 1, x^m - 1, y^p - 1, z^r - 1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x + y);

julia> B = S(1 + y + z);

julia> C = S(1 + w + z);

julia> D = S(1 + w + x);

julia> c = TetravariateTetracycleCode(l, m, p, r, A, B, C, D);

julia> code_n(c), code_k(c)
(486, 12)
```

Here is the novel `[[648, 18, d]]` tetravariate tetracycle code.

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p, r = 4, 3, 3, 3;

julia> R, (w, x, y, z) = polynomial_ring(GF(2), [:w, :x, :y, :z]);

julia> I = ideal(R, [w^l - 1, x^m - 1, y^p - 1, z^r - 1]);

julia> S, _ = quo(R, I);

julia> A = S((1 + x^2 )*(1 + w*x*y*z^2));

julia> B = S((1 + x^2)*(1 + w*x^3*y^2*z));

julia> C = S(1 + w^2*x^2*y^2*z^2);

julia> D = S(1 + w^3*x^3*y^3*z^3);

julia> c = TetravariateTetracycleCode(l, m, p, r, A, B, C, D);

julia> code_n(c), code_k(c)
(648, 18)
```

### Fields
    $TYPEDFIELDS
"""
struct TetravariateTetracycleCode <: AbstractCSSCode
     """Order of the first abelian group in ``\\mathbb{F}_2[\\mathbb{Z}_\\ell \\times \\mathbb{Z}_m \\times \\mathbb{Z}_p \\times \\mathbb{Z}_r]``"""
    l::Int
    """Order of the second abelian group in ``\\mathbb{F}_2[\\mathbb{Z}_\\ell \\times \\mathbb{Z}_m \\times \\mathbb{Z}_p \\times \\mathbb{Z}_r]``"""
    m::Int
    """Order of the third abelian group in ``\\mathbb{F}_2[\\mathbb{Z}_\\ell \\times \\mathbb{Z}_m \\times \\mathbb{Z}_p \\times \\mathbb{Z}_r]``"""
    p::Int
    """Order of the fourth abelian group in ``\\mathbb{F}_2[\\mathbb{Z}_\\ell \\times \\mathbb{Z}_m \\times \\mathbb{Z}_p \\times \\mathbb{Z}_r]``"""
    r::Int
     """First tetravariate polynomial in quotient ring ``\\frac{\\mathbb{F}_2[w, x, y, z]}{\\langle w^\\ell-1, x^m-1, y^p-1 \\rangle, y^p-1 \\rangle}``"""
    A::MPolyQuoRingElem{FqMPolyRingElem}
         """First tetravariate polynomial in quotient ring ``\\frac{\\mathbb{F}_2[w, x, y, z]}{\\langle w^\\ell-1, x^m-1, y^p-1 \\rangle, y^p-1 \\rangle}``"""
    B::MPolyQuoRingElem{FqMPolyRingElem}
         """First tetravariate polynomial in quotient ring ``\\frac{\\mathbb{F}_2[w, x, y, z]}{\\langle w^\\ell-1, x^m-1, y^p-1 \\rangle, y^p-1 \\rangle}``"""
    C::MPolyQuoRingElem{FqMPolyRingElem}
         """First tetravariate polynomial in quotient ring ``\\frac{\\mathbb{F}_2[w, x, y, z]}{\\langle w^\\ell-1, x^m-1, y^p-1 \\rangle, y^p-1 \\rangle}``"""
    D::MPolyQuoRingElem{FqMPolyRingElem}

    function TetravariateTetracycleCode(l, m, p, r, A, B, C, D)
        l > 0 || throw(ArgumentError("l must be positive"))
        m > 0 || throw(ArgumentError("m must be positive"))
        p > 0 || throw(ArgumentError("p must be positive"))
        r > 0 || throw(ArgumentError("r must be positive"))
        Rₒ = parent(A)
        R, (w, x, y, z) = polynomial_ring(GF(2), [:w, :x, :y, :z])
        I = ideal(R, [w^l-1, x^m-1, y^p-1, z^r-1])
        Rₑₓₚ, _ = quo(R, I)
        if base_ring(Rₒ) != base_ring(Rₑₓₚ) || modulus(Rₒ) != modulus(Rₑₓₚ)
            throw(ArgumentError("Polynomials must be in R/⟨w^$l-1, x^$m-1, y^$p-1, z^$r-1⟩"))
        end
        parent(B) != Rₒ && throw(ArgumentError("All polynomials must be in the same quotient ring"))
        parent(C) != Rₒ && throw(ArgumentError("All polynomials must be in the same quotient ring"))
        parent(D) != Rₒ && throw(ArgumentError("All polynomials must be in the same quotient ring"))
        new(l, m, p, r, A, B, C, D)
    end
end

function _gf2_to_int(M_gf2)
    m, n = size(M_gf2)
    M_int = zeros(Int, m, n)
    for i in 1:m, j in 1:n
        M_int[i, j] = iszero(M_gf2[i, j]) ? 0 : 1
    end
    return M_int
end

function _polynomial_to_circulant_matrix(f, l, m, p, r)
    f_lift = lift(f)
    n = l*m*p*r
    M = zero_matrix(GF(2), n, n)
    for i in 0:l-1, j in 0:m-1, k in 0:p-1, q in 0:r-1
        col_idx = i*(m*p*r)+j*(p*r)+k*r+q+1
        for term in terms(f_lift)
            c  = coeff(term, 1)
            i₂ = degree(term, 1)
            j₂ = degree(term, 2)
            k₂ = degree(term, 3)
            q₂ = degree(term, 4)
            i₃ = mod(i₂+i, l)
            j₃ = mod(j₂+j, m)
            k₃ = mod(k₂+k, p)
            q₃ = mod(q₂+q, r)
            row_idx = i₃*(m*p*r)+j₃*(p*r)+k₃*r+q₃+1
            M[row_idx, col_idx] += c
        end
    end
    return M
end

function boundary_maps(c::TetravariateTetracycleCode)
    l, m, p, r = c.l, c.m, c.p, c.r
    M_A = _polynomial_to_circulant_matrix(c.A, l, m, p, r)
    M_B = _polynomial_to_circulant_matrix(c.B, l, m, p, r)
    M_C = _polynomial_to_circulant_matrix(c.C, l, m, p, r)
    M_D = _polynomial_to_circulant_matrix(c.D, l, m, p, r)
    n = size(M_A, 1)
    zero_block = zero_matrix(GF(2), n, n)
    ∂₄ = vcat(M_A, M_B, M_C, M_D)
    ∂₃ = hvcat((4, 4, 4, 4, 4, 4),
        M_B        , M_A        , zero_block, zero_block,
        M_C        , zero_block , M_A       , zero_block,
        M_D        , zero_block , zero_block, M_A       ,
        zero_block , M_C        , M_B       , zero_block,
        zero_block , M_D        , zero_block, M_B       ,
        zero_block , zero_block , M_D       , M_C
    )
    @assert iszero(∂₃*∂₄)
    ∂₂ = hvcat((6, 6, 6, 6),
        zero_block, zero_block, zero_block, M_D       , M_C       , M_B       ,
        zero_block, M_D       , M_C       , zero_block, zero_block, M_A       ,
        M_D       , zero_block, M_B       , zero_block, M_A       , zero_block,
        M_C       , M_B       , zero_block, M_A       , zero_block, zero_block
    )
    @assert iszero(∂₂*∂₃)
    ∂₁ = hcat(M_A, M_B, M_C, M_D)
    @assert iszero(∂₁*∂₂)
    return _gf2_to_int(∂₁), _gf2_to_int(∂₂), _gf2_to_int(∂₃), _gf2_to_int(∂₄)
end

function parity_matrix_xz(c::TetravariateTetracycleCode)
    _, ∂₂, ∂₃, _ = boundary_maps(c)
    return ∂₂, Matrix(∂₃')
end

parity_matrix_x(c::TetravariateTetracycleCode) = boundary_maps(c)[2]

parity_matrix_z(c::TetravariateTetracycleCode) = Matrix(boundary_maps(c)[3]')

metacheck_matrix_z(c::TetravariateTetracycleCode) = Matrix(boundary_maps(c)[4]')

metacheck_matrix_x(c::TetravariateTetracycleCode) = boundary_maps(c)[1]
