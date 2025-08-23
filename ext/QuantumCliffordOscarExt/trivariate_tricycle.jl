"""
    $TYPEDEF

A quantum CSS code constructed from three trivariate polynomials over a finite field
from [jacob2025singleshotdecodingfaulttolerantgates](@cite).

#### Example

Here is the `[[72, 6, 6]]` trivariate tricycle code from Table I from  [jacob2025singleshotdecodingfaulttolerantgates](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p = 4, 3, 2;

julia> R, (x, y, z) = polynomial_ring(GF(2), [:x, :y, :z]);

julia> I = ideal(R, [x^l - 1, y^m - 1, z^p - 1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + y + x*y^2);

julia> B = S(1 + y*z + x^2*y^2);

julia> C = S(1 + x*y^2*z + x^2*y);

julia> c = TrivariateTricycleCode(l, m, p, A, B, C);

julia> code_n(c), code_k(c)
(72, 6)
```

Here is the `[[432, 12, 12]` trivariate tricycle code from Table I from  [jacob2025singleshotdecodingfaulttolerantgates](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p = 6, 6, 4;

julia> R, (x, y, z) = polynomial_ring(GF(2), [:x, :y, :z]);

julia> I = ideal(R, [x^l - 1, y^m - 1, z^p - 1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x*y*z^3 + x^3*y^4*z^2);

julia> B = S(1 + x^3*y*z^2 + x^3*y^2*z^3);

julia> C = S(1 + x^4*y^3*z^3 + x^5*z^2);

julia> c = TrivariateTricycleCode(l, m, p, A, B, C);

julia> code_n(c), code_k(c)
(432, 12)
```

### Fields
    $TYPEDFIELDS
"""
struct TrivariateTricycleCode <: AbstractCSSCode
    """Size of the first cyclic dimension"""
    l::Int
    """Size of the second cyclic dimension"""
    m::Int
    """Size of the third cyclic dimension"""
    p::Int
    """First trivariate polynomial in quotient ring R/⟨xˡ-1, yᵐ-1, zᵖ-1⟩"""
    A::MPolyQuoRingElem{FqMPolyRingElem}
    """Second trivariate polynomial in quotient ring R/⟨xˡ-1, yᵐ-1, zᵖ-1⟩"""
    B::MPolyQuoRingElem{FqMPolyRingElem}
    """Third trivariate polynomial in quotient ring R/⟨xˡ-1, yᵐ-1, zᵖ-1⟩"""
    C::MPolyQuoRingElem{FqMPolyRingElem}

    function TrivariateTricycleCode(l::Int, m::Int, p::Int, A::MPolyQuoRingElem{FqMPolyRingElem}, B::MPolyQuoRingElem{FqMPolyRingElem}, C::MPolyQuoRingElem{FqMPolyRingElem})
        l > 0 || throw(ArgumentError("l must be positive"))
        m > 0 || throw(ArgumentError("m must be positive"))
        p > 0 || throw(ArgumentError("p must be positive"))
        Rₒ = parent(A)
        R, (x,y,z) = polynomial_ring(GF(2), [:x, :y, :z])
        I = ideal(R, [x^l-1, y^m-1, z^p-1])
        Rₑₓₚ, _ = quo(R, I)
        base_ring(Rₒ) != base_ring(Rₑₓₚ) && throw(ArgumentError("A must be in R/⟨x^$l-1, y^$m-1, z^$p-1⟩"))
        modulus(Rₒ) != modulus(Rₑₓₚ) && throw(ArgumentError("A must be in R/⟨x^$l-1, y^$m-1, z^$p-1⟩"))
        parent(B) != Rₒ && throw(ArgumentError("B must be in same ring as A"))
        parent(C) != Rₒ && throw(ArgumentError("C must be in same ring as A"))
        new(l, m, p, A, B, C)
    end
end

"""Convert a matrix over GF(2) to an integer matrix."""
function _gf2_to_int_mat(M_gf2)
    m, n = size(M_gf2)
    Mᵢₙₜ = zeros(Int, m, n)
    for i in 1:m, j in 1:n
        elem = M_gf2[i, j]
        Mᵢₙₜ[i, j] = iszero(elem) ? 0 : 1
    end
    return Mᵢₙₜ
end

"""Construct a 3D circulant matrix from a trivariate polynomial."""
function _circulant_matrix_3d(poly, l, m, p)
    n = l*m*p
    M = zero_matrix(GF(2), n, n)
    for i in 0:l-1, j in 0:m-1, k in 0:p-1
        col_idx = i*(m*p)+j*p+k+1
        for term in terms(poly)
            c = coeff(term, 1)
            i₂ = degree(term, 1)
            j₂ = degree(term, 2)
            k₂ = degree(term, 3)
            i₃ = mod(i₂+i, l)
            j₃ = mod(j₂+j, m)
            k₃ = mod(k₂+k, p)
            row_idx = i₃*(m*p)+j₃*p+k₃+1
            M[row_idx, col_idx] += c
        end
    end
    return M
end

"""Convert a quotient ring polynomial to its circulant matrix representation."""
function _polynomial_to_circulant_matrix(f, l, m, p)
    f_lift = lift(f)
    return _circulant_matrix_3d(f_lift, l, m, p)
end

"""Compute the transpose of a polynomial in the group algebra. For a
polynomial f = ∑ cᵢⱼₖ xⁱyʲzᵏ, the transpose is fᵀ = ∑ cᵢⱼₖ x⁻ⁱy⁻ʲz⁻ᵏ."""
function _polynomial_transpose(f, l, m, p)
    f_lift = lift(f)
    Rₒ = parent(f_lift)
    x, y, z = gens(Rₒ)
    fₜᵣₐₙₛ = zero(Rₒ)
    monoms = monomials(f_lift)
    coeffs = coefficients(f_lift)
    for (mono, c) in zip(monoms, coeffs)
        i = degree(mono, 1)
        j = degree(mono, 2)
        k = degree(mono, 3)
        iᵢₙᵥ = mod(-i, l)
        jᵢₙᵥ = mod(-j, m)
        kᵢₙᵥ = mod(-k, p)
        monoᵢₙᵥ = x^iᵢₙᵥ*y^jᵢₙᵥ*z^kᵢₙᵥ
        fₜᵣₐₙₛ += c*monoᵢₙᵥ
    end
    return Rₒ(fₜᵣₐₙₛ)
end

function boundary_maps(c::TrivariateTricycleCode)
    l, m, p = c.l, c.m, c.p
    A, B, C = c.A, c.B, c.C
    M_A = _polynomial_to_circulant_matrix(A, l, m, p)
    M_B = _polynomial_to_circulant_matrix(B, l, m, p)
    M_C = _polynomial_to_circulant_matrix(C, l, m, p)
    Aₜᵣₐₙₛ = _polynomial_transpose(A, l, m, p)
    Bₜᵣₐₙₛ = _polynomial_transpose(B, l, m, p)
    Cₜᵣₐₙₛ = _polynomial_transpose(C, l, m, p)
    M_Aₜᵣₐₙₛ = _polynomial_to_circulant_matrix(Aₜᵣₐₙₛ, l, m, p)
    M_Bₜᵣₐₙₛ = _polynomial_to_circulant_matrix(Bₜᵣₐₙₛ, l, m, p)
    M_Cₜᵣₐₙₛ = _polynomial_to_circulant_matrix(Cₜᵣₐₙₛ, l, m, p)
    n_block = l*m*p
    zero_block = zero_matrix(GF(2), n_block, n_block)
    H_X = hcat(M_A, M_B, M_C) # Eq. 14
    H_Z = vcat(
        hcat(zero_block, M_Cₜᵣₐₙₛ, M_Bₜᵣₐₙₛ),
        hcat(M_Cₜᵣₐₙₛ, zero_block, M_Aₜᵣₐₙₛ),
        hcat(M_Bₜᵣₐₙₛ, M_Aₜᵣₐₙₛ, zero_block)
    ) # Eq. 13
    M_Z = hcat(M_Aₜᵣₐₙₛ, M_Bₜᵣₐₙₛ, M_Cₜᵣₐₙₛ) # Eq. 12
    H_X = _gf2_to_int_mat(H_X)
    H_Z = _gf2_to_int_mat(H_Z)
    M_Z = _gf2_to_int_mat(M_Z)
    return H_X, H_Z, M_Z
end

parity_matrix_x(c::TrivariateTricycleCode) = boundary_maps(c)[1]

parity_matrix_z(c::TrivariateTricycleCode) = boundary_maps(c)[2]

metacheck_matrix_z(c::TrivariateTricycleCode) = boundary_maps(c)[3]
