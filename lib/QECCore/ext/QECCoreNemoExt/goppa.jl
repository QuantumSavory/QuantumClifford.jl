"""The family of Goppa codes, as discovered by Denisovich Goppa, in his 1970 paper
[goppa1970new](@cite). An irreducible binary Goppa code is characterized by two key elements:

- Goppa Polynomial `(g(x))`: A polynomial of degree `t` defined over the finite field ``\\mathbb{F}_{2^m}``. This polynomial must have no repeated roots. In other words, Let g(x) be a monic polynomial over the Fₓᵐ and let L={γ₀, ..., γₙ₋₁} be a set of n elements of Fₓᵐ such that g(γᵢ) ≠ 0 for 0 ≤ i < n.
- Support List `(L)`: A list containing `n` distinct elements from the finite field ``\\mathbb{F}_{2^m}``. These elements must not be roots of the Goppa polynomial g(x). 

The sequence L plays a crucial role in defining the structure of an irreducible binary Goppa code. The steps involved in forming L are as follows:

- Field Selection: A finite field, ``\\mathbb{F}_{2^m}``, is chosen based on the extension degree `m`.
- Goppa Polynomial: The input polynomial g(x) is used, which should be irreducible and have no repeated roots within ``\\mathbb{F}_{2^m}``.
- Elements Selection for `L`: The sequence `L` consists of `n` distinct elements chosen from ``\\mathbb{F}_{2^m}``. Importantly, no element in L can be a root of g(x).

The parity-check matrix `(H)` of an irreducible binary Goppa code can be broken down into product of three simpler matrices:

- Vandermonde Matrix `(V)`:
  - Each row of `V` is constructed using the elements from the support list `(L)` raised to different powers.

```math
\\begin{matrix}
1                & 1                  & 1                & \\cdots & 1 \\\\
(L_1)^1          & (L_2)^1            & (L_3)^1          & \\cdots & (L_n)^1 \\\\
(L_1)^2          & (L_2)^2            & (L_3)^2          & \\cdots & (L_n)^2 \\\\
\\vdots          & \\vdots            & \\vdots          & \\ddots & \\vdots \\\\
(L_1)^{t-1}      & (L_2)^{t-1}        & (L_3)^{t-1}      & \\cdots & (L_n)^{t-1}
\\end{matrix}
```

- Diagonal Matrix `(D)`:
  - Each diagonal element in `D` is related to the roots of the Goppa polynomial `(g(x))`.

```math
\\begin{bmatrix}
\\frac{1}{g(L_1)} & 0               & 0                 & \\cdots & 0 \\\\
0               & \\frac{1}{g(L_2)} & 0                 & \\cdots & 0 \\\\
0               & 0                 & \\frac{1}{g(L_3)} & \\cdots & 0 \\\\
\\vdots         & \\vdots           & \\vdots           & \\ddots & \\vdots \\\\
0               & 0                 & 0                 & \\cdots & \\frac{1}{g(L_n)}
\\end{bmatrix}
```

- X Matrix `(X)`:
  - Each element in `X` is related to the coefficients of the Goppa polynomial `(g(x))`

```math
\\begin{bmatrix}
g_t           & 0             & 0             & \\cdots & 0 \\\\
g_{t-1}       & g_t           & 0             & \\cdots & 0 \\\\
g_{t-2}       & g_{t-1}       & g_t           & \\cdots & 0 \\\\
\\vdots       & \\vdots       & \\vdots       & \\ddots & \\vdots \\\\
g_1           & g_2           & g_3           & \\cdots & g_t
\\end{bmatrix}
```

The matrix `VD` can be viewed as parity check matrix for `Γ(L, g)` over ``\\mathbb{F}_{2^m}``. The matrix is as follows:

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
    ``\\mathbb{F}_{2^m}`` can be represented by a `m`-tuple/binary column vector of length `m`.
    If each entry of `H` is replaced by its corresponding `m`-tuple/binary column vector of length
    `m` over ``\\mathbb{F}_{2^m}``, we obtain a binary parity check matrix for the goppa code.

You might be interested in consulting [berlekamp1973goppa](@cite), [mceliece1978public](@cite),
[patterson1975algebraic](@cite), [sugiyama1975method](@cite), [van1988classical](@cite),
[bernstein2008attacking](@cite), [wirtz1988parameters](@cite) and [singh2019code](@cite)
an as well. 

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/gappa).

!!! note
    In Nemo, taking a random monic poly of degree `n`, this poly is irreducible with
    probability `1/n`. One in `n` monic polynomials is, on average irreducible. To increase
    probability of success of getting irreducible polynomial, use more iterations.
"""
abstract type AbstractPolynomialCode <: AbstractCECC end

struct GoppaCode <: AbstractPolynomialCode
    m::Int 
    t::Int 
    seed::Int

    function GoppaCode(m, t; seed::Int=42)
        if m < 3 || t < 2 || t >= 2 ^ (m - 1)
            throw(ArgumentError("Invalid parameters: 'm' and 't' must be positive. Additionally, 'm' is >= 3 and t >= 2 to obtain a valid code and to tractable."))
        end
        new(m, t, seed)
    end
end

function generator_polynomial(ga::GoppaCode)
    GF2ʳ, o = finite_field(2, ga.m, :o)
    k = GF(2, ga.m)
    po, b = polynomial_ring(k)
    gx = FqPolyRingElem
    rng = MersenneTwister(ga.seed)
    for i in 1:500
        gx = rand(rng, po, 1:ga.t) + b^ga.t
        if is_irreducible(gx) == true
            return gx  
        end
    end
end

function parity_matrix(ga::GoppaCode)
    n = 2^ga.m
    F, o = finite_field(2, ga.m, :o)
    g = generator_polynomial(ga)
    # L = {oⁱ | g(oⁱ) ≠ 0}
    L = [o^i for i in 0:n-1 if evaluate(g, o^i) != 0]
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
