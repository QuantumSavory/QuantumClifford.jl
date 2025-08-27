"""
$TYPEDEF

Generalized bicycle codes ([koukoulekidis2024smallquantumcodesalgebraic](@cite))

A generalized bicycle quantum LDPC code constructed from two polynomials in ``\\mathbb{F}_2[x]/(x^l - 1)``.

Here is an example of a [[10, 2, 3]] GB code from the Appendix B of
[koukoulekidis2024smallquantumcodesalgebraic](@cite) with lift size of 5
build out of two polymonials in ``\\mathbb{F}_2[x]/(x^l - 1)``.

```jldoctest
julia> import Hecke: polynomial_ring, GF; using QuantumClifford.ECC;

julia> R, x = polynomial_ring(GF(2), "x");

julia> l = 5;

julia> a = 1 + x^4;

julia> b = 1 + x + x^2 + x^4;

julia> c = GeneralizedBicycleCode(a, b, l);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(10, 2, 3)
```

Here is an example of a [[12, 2, 3]] GB code from the Appendix B of
[koukoulekidis2024smallquantumcodesalgebraic](@cite) with lift size of 6
build out of two polymonials in ``\\mathbb{F}_2[x]/(x^l - 1)``.

```jldoctest
julia> import Hecke: polynomial_ring, GF; using QuantumClifford.ECC;

julia> R, x = polynomial_ring(GF(2), "x");

julia> l = 6;

julia> a = 1 + x + x^2 + x^5;

julia> b = 1 + x + x^3 + x^5;

julia> c = GeneralizedBicycleCode(a, b, l);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(12, 2, 3)
```

See also: [`generalized_bicycle_codes`](@ref)

### Fields
    $TYPEDFIELDS
"""
struct GeneralizedBicycleCode <: AbstractCSSCode
    """First generator polynomial in 𝔽₂[x]/(xˡ - 1)."""
    a::FqPolyRingElem
    """Second generator polynomial in 𝔽₂[x]/(xˡ - 1)."""
    b::FqPolyRingElem
    """The lift size which corresponds to dimension of circulant matrices."""
    l::Int
    function GeneralizedBicycleCode(a::FqPolyRingElem, b::FqPolyRingElem, l::Int)
        l <= 0 && throw(ArgumentError("Block length must be positive."))
        (base_ring(a) != base_ring(b)) && throw(ArgumentError("Polynomials must be from the same ring."))
        (characteristic(base_ring(a)) != 2) && throw(ArgumentError("Polynomials must be over 𝔽₂"))   
        (degree(a) >= l || degree(b) >= l) && throw(ArgumentError("Polynomial degrees must be < l."))
        new(a, b, l)
    end
end

function parity_matrix_xz(c::GeneralizedBicycleCode)
    A = circulant_matrix_from_polynomial_ring(c.l, c.a)
    B = circulant_matrix_from_polynomial_ring(c.l, c.b)
    hx = hcat(A, B)
    hz = hcat(B', A')
    return hx, hz
end

parity_matrix_x(c::GeneralizedBicycleCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::GeneralizedBicycleCode) = parity_matrix_xz(c)[2]

code_n(c::GeneralizedBicycleCode) = 2*c.l

code_k(c::GeneralizedBicycleCode) = 2*degree(gcd(c.a, c.b, gen(parent(c.a))^c.l - 1))
