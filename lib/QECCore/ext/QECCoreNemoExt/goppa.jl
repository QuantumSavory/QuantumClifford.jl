abstract type AbstractPolynomialCode <: AbstractCECC end

"""
    $TYPEDEF

The family of Goppa codes, as discovered by Denisovich Goppa, in his 1970 paper [goppa1970new](@cite). The
binary Goppa code is characterized by two key elements:

- Goppa polynomial: A *monic* polynomial of degree t over ``\\mathbb{F}_{2^m}`` with no repeated roots. For a fixed support set ``L = {\\gamma_0, \\dots, \\gamma_{n-1}} \\subseteq \\mathbb{F}_{2^m}``, the polynomial satisfies ``g(\\gamma_i) \\neq 0`` for all ``0 \\leq i < n``.  
- Support list: A list of n distinct elements ``{\\gamma_0, \\dots, \\gamma_{n-1}}`` from ``\\mathbb{F}_{2^m}`` such that ``g(\\gamma_i) \\neq 0`` for all ``0 \\leq i < n`` (i.e., none are roots of the Goppa polynomial g(x).

The set ``L = {\\gamma_0, \\dots, \\gamma_{n-1}}`` defines the code’s structure, with n distinct
elements from ``\\mathbb{F}_{2^m}``. The Goppa polynomial g(x) is irreducible over ``\\mathbb{F}_{2^m}``
and satisfy ``g(\\gamma_i) \\neq 0)`` for all ``(\\gamma_i \\in L``. The parity-check matrix H of
an binary Goppa code can be broken down into product of three simpler matrices:

## Vandermonde Matrix

Each row of V is constructed using the elements from the support list L raised to different powers.

```math
\\begin{matrix}
1                & 1                  & 1                & \\cdots & 1 \\\\
(L_1)^1          & (L_2)^1            & (L_3)^1          & \\cdots & (L_n)^1 \\\\
(L_1)^2          & (L_2)^2            & (L_3)^2          & \\cdots & (L_n)^2 \\\\
\\vdots          & \\vdots            & \\vdots          & \\ddots & \\vdots \\\\
(L_1)^{t-1}      & (L_2)^{t-1}        & (L_3)^{t-1}      & \\cdots & (L_n)^{t-1}
\\end{matrix}
```

## Diagonal Matrix

Each diagonal element in D is related to the roots of the Goppa polynomial g(x).

```math
\\begin{bmatrix}
\\frac{1}{g(L_1)} & 0               & 0                 & \\cdots & 0 \\\\
0               & \\frac{1}{g(L_2)} & 0                 & \\cdots & 0 \\\\
0               & 0                 & \\frac{1}{g(L_3)} & \\cdots & 0 \\\\
\\vdots         & \\vdots           & \\vdots           & \\ddots & \\vdots \\\\
0               & 0                 & 0                 & \\cdots & \\frac{1}{g(L_n)}
\\end{bmatrix}
```

## X Matrix

Each element in X is related to the coefficients of the Goppa polynomial g(x).

```math
\\begin{bmatrix}
g_t           & 0             & 0             & \\cdots & 0 \\\\
g_{t-1}       & g_t           & 0             & \\cdots & 0 \\\\
g_{t-2}       & g_{t-1}       & g_t           & \\cdots & 0 \\\\
\\vdots       & \\vdots       & \\vdots       & \\ddots & \\vdots \\\\
g_1           & g_2           & g_3           & \\cdots & g_t
\\end{bmatrix}
```

The matrix V*D can be viewed as parity check matrix for Γ(L, g) over ``\\mathbb{F}_{2^m}``:

```math
\\begin{bmatrix}
1 \\cdot g(L_1)^{-1}       & 1 \\cdot g(L_2)^{-1}       & 1 \\cdot g(L_3)^{-1}       & \\cdots & 1 \\cdot g(L_n)^{-1} \\\\
L_1 \\cdot g(L_1)^{-1}     & L_2 \\cdot g(L_2)^{-1}     & L_3 \\cdot g(L_3)^{-1}     & \\cdots & L_n \\cdot g(L_n)^{-1} \\\\
L_1^2 \\cdot g(L_1)^{-1}   & L_2^2 \\cdot g(L_2)^{-1}   & L_3^2 \\cdot g(L_3)^{-1}   & \\cdots & L_n^2 \\cdot g(L_n)^{-1} \\\\
\\vdots                    & \\vdots                    & \\vdots                    & \\ddots & \vdots \\\\
L_1^{t-1} \\cdot g(L_1)^{-1} & L_2^{t-1} \\cdot g(L_2)^{-1} & L_3^{t-1} \\cdot g(L_3)^{-1} & \\cdots & L_n^{t-1} \\cdot g(L_n)^{-1}
\\end{bmatrix}
```

!!! note
    The entries of matrix over field elements are in ``\\mathbb{F}_{2^m}``. Each element in
    ``\\mathbb{F}_{2^m}`` can be represented by a m-tuple/binary column vector of length m.
    If each entry of H is replaced by its corresponding m-tuple/binary column vector of length
    m over ``\\mathbb{F}_{2^m}``, we obtain a binary parity check matrix for the goppa code.

You might be interested in consulting [berlekamp1973goppa](@cite), [mceliece1978public](@cite),
[patterson1975algebraic](@cite), [sugiyama1975method](@cite), [van1988classical](@cite),
[bernstein2008attacking](@cite), [wirtz1988parameters](@cite) and [singh2019code](@cite)
an as well. 

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/gappa).

!!! note
    In Nemo, taking a random monic poly of degree `n`, this poly is irreducible with
    probability `1/n`. One in `n` monic polynomials is, on average irreducible. To increase
    probability of success of getting irreducible polynomial, use more iterations.

### Example

Here is an example of `[8, 2, d]` Goppa code from [slides](https://crypto-kantiana.com/elena.kirshanova/talks/Talk_McEliece.pdf).

```jldoctest
julia> using QECCore; using Nemo

julia> t = 2;

julia> m = 3;

julia> F, α = finite_field(2, m, :α);

julia> R, x = polynomial_ring(F, :x);

julia> g = x^t + x + 1;

julia> L = [F(0), F(1), α, α^2, α + 1, α^2 + α, α^2 + α + 1, α^2 + 1];

julia> c = GoppaCode(m, t, g, L);

julia> code_n(c), code_k(c)
(8, 2)
```

If no support set L is specified, it defaults to all non-roots of g(x) in ``\\mathbb{F}_{2^m}``

```math
\\begin{aligned}
L = { \\alpha \\in \\mathbb{F}_{2^m} \\mid g(\\alpha) \\neq 0}
\\end{aligned}
```

```jldoctest
julia> using QECCore; using Nemo

julia> m, t =  5, 4;

julia> c = GoppaCode(m, t);

julia> c.g
x^4 + α^4*x^3 + (α^3 + α^2)*x + α^3 + α

julia> code_n(c), code_k(c)
(32, 12)
```

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/goppa).

### Fields
    $TYPEDFIELDS
"""
struct GoppaCode <: AbstractPolynomialCode
    """ The extension degree of ``\\mathbb{F}_{2^m}`` which determines the size of the field
    (2^m elements) and the binary expansion length."""
    m::Int
    """The degree of the Goppa polynomial `g(x)` which controls the error-correction capability."""
    t::Int
    """The Goppa polynomial over ``\\mathbb{F}_{2^m}`` which must satisfy: degree(g) = t, and
    g(αᵢ) ≠ 0 ∀ αᵢ ∈ L."""
    g::FqPolyRingElem
    """The support set of distinct elements from ``\\mathbb{F}_{2^m}`` which defines code length
    n = length(L). It must satisfy g(α) ≠ 0 ∀ α ∈ L."""
    L::Vector{FqFieldElem}
    """Random number generator seed for reproducible polynomial generation."""
    seed::Int

    function GoppaCode(m, t, g::FqPolyRingElem, L::Vector{FqFieldElem}; seed::Int=42)
        (m < 3 || t < 2 || t >= 2^(m - 1)) && throw(ArgumentError("m ≥ 3 and t ≥ 2 required, with t < 2^(m-1)"))
        degree(g) != t && throw(ArgumentError("The Goppa polynomial must have degree t"))
        !is_monic(g) && throw(ArgumentError("The Goppa polynomial must be monic."))
        any(evaluate(g, α) == 0 for α in L) && throw(ArgumentError("L cannot contain roots of the Goppa polynomial."))
        length(L) - m*t ≤ 0 && throw(ArgumentError("Parameters would produce invalid code: length(L) - m*t = $(length(L)) - $m*$t ≤ 0."))
        new(m, t, g, L, seed)
    end
end
    
function GoppaCode(m, t, g::FqPolyRingElem; seed::Int=42)
    F, α = finite_field(2, m, :α)
    L = [a for a in F if evaluate(g, a) != 0]
    return GoppaCode(m, t, g, L; seed)
end

function GoppaCode(m, t; seed=42)
    F, α = finite_field(2, m, :α)
    R, x = polynomial_ring(F, "x")
    rng = MersenneTwister(seed)
    for _ in 1:500
        coeffs = rand(rng, F, t)
        g = R([coeffs...; one(F)])
        if is_irreducible(g)
            L = [a for a in F if evaluate(g, a) != 0]
            return GoppaCode(m, t, g, L; seed)
        end
    end
    throw(ErrorException("Failed to find irreducible polynomial after 500 attempts"))
end

function parity_matrix(ga::GoppaCode)
    m, t, g, L = ga.m, ga.t, ga.g, ga.L
    F = parent(L[1])
    # D = diag(1/g(oᵢ))
    Dᵢₙᵥ = [inv(evaluate(g, o_i)) for o_i in L]
    # V = [oᵢʲ⁻¹]
    V = [o_i^(j-1) for j in 1:ga.t, o_i in L]
    H_poly = [V[i,j] * Dᵢₙᵥ[j] for i in 1:ga.t, j in 1:length(L)]
    H = Matrix{Bool}(undef, ga.m * ga.t, length(L))
    for i in 1:ga.t
        row_range = (i-1)*ga.m+1 : i*ga.m
        for j in 1:length(L)
            H[row_range, j] .= [coeff(H_poly[i,j], k) != 0 for k in 0:ga.m-1]
        end
    end
    return H
end

generator_polynomial(c::GoppaCode) = c.g
