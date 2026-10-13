"""
The family of Bose–Chaudhuri–Hocquenghem (BCH) codes, as discovered in 1959 by
Alexis Hocquenghem [hocquenghem1959codes](@cite), and independently in 1960 by
Raj Chandra Bose and D.K. Ray-Chaudhuri [bose1960class](@cite).

The BCH code, denoted as `BCH(m, t)`, is defined by `m`, a positive integer that
specifies the degree of the finite (Galois) field `GF(2ᵐ)`, and `t`, a positive
integer that indicates the number of correctable errors. The binary parity check
matrix can be obtained from the following matrix over ``GF(2ᵐ)`` field elements:

```math
\\begin{aligned}
\\begin{matrix}
1 & (\\alpha^1)^1 & (\\alpha^1)^2 & (\\alpha^1)^3 & \\dots & (\\alpha^1)^{n-1} \\\\
1 & (\\alpha^3)^1 & (\\alpha^3)^2 & (\\alpha^3)^3 & \\dots & (\\alpha^3)^{n-1} \\\\
1 & (\\alpha^5)^1 & (\\alpha^5)^2 & (\\alpha^5)^3 & \\dots & (\\alpha^5)^{n-1} \\\\
\\vdots & \\vdots & \\vdots & \\vdots & \\ddots & \\vdots \\\\
1 & (\\alpha^{2t-1})^1 & (\\alpha^{2t-1})^2 & (\\alpha^{2t-1})^3 & \\dots & (\\alpha^{2t-1})^{n-1}
\\end{matrix}
\\end{aligned}
```

The entries of the matrix are in `GF(2ᵐ)`. Each element in `GF(2ᵐ)` can be represented
by an `m`-tuple (a binary column vector of length `m`). If each entry of `H` is replaced
by its corresponding `m`-tuple, we obtain a binary parity check matrix for the code.

The BCH code is cyclic as its generator polynomial, `g(x)` divides `xⁿ - 1`, so
`mod (xⁿ - 1, g(x)) = 0`.

You might be interested in consulting [bose1960further](@cite) and [error2024lin](@cite)
as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/q-ary_bch).
"""
struct BCH <: AbstractPolynomialCode
    m::Int
    t::Int
    function BCH(m, t)
        m < 3 && throw(ArgumentError("m must be greater than or equal to 3"))
        t >= 2^(m - 1) && throw(ArgumentError("t must be less than 2ᵐ ⁻ ¹"))
        m * t > 2^m - 1 && throw(ArgumentError("m*t must be less than or equal to 2ᵐ - 1"))
        new(m, t)
    end
end

code_n(b::BCH) = 2 ^ b.m - 1
