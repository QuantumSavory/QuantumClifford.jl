"""
The extended generalized bicycle code is a family of quantum LDPC codes generated
through *algebraic extension* of a base [`GeneralizedBicycleCode`](@ref). Starting
with initial generating polynomials ``a(x), b(x) \\in \\mathbb{F}_2^{\\langle\\ell\\rangle}``,
the extended GB codes are constructed by polynomial multiplication, where for each extension
step ``m``, an extension polynomial 

```math
\\begin{aligned}
p^{(m)}(x) \\in \\mathbb{F}_2^{\\langle\\kappa_m\\ell\\rangle}
\\end{aligned}
```

is selected to produce extended polynomials:

```math
\\begin{aligned}
a^{(m)}(x) &= p^{(m)}(x)a(x), \\\\
b^{(m)}(x) &= p^{(m)}(x)b(x)
\\end{aligned}
```

These extended polynomials are then used to form ``\\kappa_m\\ell \\times \\kappa_m\\ell`` circulant
matrices ``A_m`` and ``B_m``, which are combined into a parity-check matrix ``H_m`` with block structure:

```math
\\begin{aligned}
H_m = \\begin{pmatrix} 
A_m|B_m & 0 \\\\ 
0 & B_m^\\top|A_m^\\top 
\\end{aligned}
```

while maintaining dimensions:

```math
\\begin{aligned}
k_m = 2\\deg\\left(\\gcd\\left(a^{(m)}, b^{(m)}, x^{\\kappa_m\\ell}-1\\right)\\right)
\\end{aligned}
```

that are always bounded below by the base code dimension `k`.

#### Example

```jldoctest
julia> import Hecke: polynomial_ring, GF, one, gen; using QuantumClifford.ECC;

julia> R, x = polynomial_ring(GF(2), "x");

julia> l = 5;

julia> a = 1 + x^4;

julia> b = 1 + x + x^2 + x^4;

julia> c = GeneralizedBicycleCode(a, b, l);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(10, 2, 3)

julia> m, p = 2, one(R);

julia> new_code = ExtendedGeneralizedBicycleCode(c, m, p);

julia> code_n(new_code), code_k(new_code), distance(new_code, DistanceMIPAlgorithm(solver=HiGHS))
(20, 2, 5)

julia> p = 1 + x;

julia> new_code = ExtendedGeneralizedBicycleCode(c, m, p);

julia> code_n(new_code), code_k(new_code), distance(new_code, DistanceMIPAlgorithm(solver=HiGHS))
(20, 4, 4)
```
"""
struct ExtendedGeneralizedBicycleCode <: AbstractCSSCode
    """The base generalized bicycle code to extend."""
    base_code::AbstractCSSCode
    """The extension index (m ≥ 1)"""
    m::Int
    """The extension polynomial ∈ 𝔽₂[((m-1)ℓ +1)]."""
    p::FqPolyRingElem
    
    function ExtendedGeneralizedBicycleCode(base_code::GeneralizedBicycleCode, m::Int, p::FqPolyRingElem)
        m ≥ 1 || throw(ArgumentError("Extension index m must be ≥ 1"))
        if m == 1
            isone(p) || throw(ArgumentError("For m=1, p must be 1"))
        else
            max_deg = (m-1)*base_code.l+1
            degree(p) < max_deg || throw(ArgumentError("For m=$m, polynomial degree must be < $max_deg"))
        end
        new(base_code, m, p)
    end
end

function parity_matrix_xz(c::ExtendedGeneralizedBicycleCode)
    ℓ = c.base_code.l
    R = parent(c.base_code.a)
    x = gen(R)
    a⁽ᵐ⁾ = mod(c.p*c.base_code.a, x^(c.m*ℓ)-1)
    b⁽ᵐ⁾ = mod(c.p*c.base_code.b, x^(c.m*ℓ)-1)
    ext_gb = GeneralizedBicycleCode(a⁽ᵐ⁾, b⁽ᵐ⁾, c.m*ℓ)
    hx, hz = parity_matrix_xz(ext_gb)
    return hx, hz
end

parity_matrix_x(c::ExtendedGeneralizedBicycleCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::ExtendedGeneralizedBicycleCode) = parity_matrix_xz(c)[2]
