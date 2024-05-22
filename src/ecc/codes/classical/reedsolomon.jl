"""The family of Reed-Solomon codes, as discovered by Reed and Solomon in their 1960 paper [reed1960polynomial](@cite). 

You might be interested in consulting [geisel1990tutorial](@cite), [wicker1999reed](@cite), [sklar2001reed](@cite), and [berlekamp1978readable](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/reed_solomon).
"""
abstract type AbstractPolynomialCode <: ClassicalCode end

struct ReedSolomon <: AbstractPolynomialCode
    m::Int
    t::Int

    function ReedSolomon(m, t)
        if m < 3 || t < 0 || t >= 2 ^ (m - 1) 
            throw(ArgumentError("Invalid parameters: m and t must be non-negative. Also, m > 3 and t < 2 ^ (m - 1) in order to obtain a valid code."))
        end
        new(m, t)
    end
end

"""
`generator_polynomial(ReedSolomon(m, t))`
- `m`: The positive integer defining the degree of the finite (Galois) field, `GF(2ᵐ)`.
- `t`: The positive integer specifying the number of correctable errors.

The generator polynomial for an RS code takes the following form:
```
g(X) = g₀ + g₁X¹ + g₂X² + ... + g₂ₜ₋₁X²ᵗ⁻¹ + X²ᵗ
```

where `X` is the indeterminate variable, `gᵢ` are the coefficients of the polynomial and `t` is the number of correctable symbol errors.

We describe the generator polynomial in terms of its `2 * t  = n - k` roots, as follows:
``` 
g(X) = (X - α¹)(X - α²)(X - α³) ... (X - α²ᵗ)
```

Degree and Parity Symbols: The degree of the generator polynomial is equal to `2 * t`, which is also the number of parity symbols added to the original data (`k` symbols) to create a codeword of length `n` `(n = k + 2 * t)`.

Roots of the Generator Polynomial: The generator polynomial has `2 * t` distinct roots, designated as `α¹, α², ... , α²ᵗ`. These roots are chosen from a Finite Galois Field. Any power of α can be used as the starting root, not necessarily `α¹` itself.

Fixed generator polynomial scheme vs variable generator polynomial scheme: Only in this construction scheme using fixed generator polynomial `g(x)`, RS codes are a subset of the Bose, Chaudhuri, and Hocquenghem (BCH) codes; hence, this relationship between the degree of the generator polynomial and the number of parity symbols holds, just as for BCH codes where degree of BCH generator polynomial, `degree(g(x)) == n - k`. Prior to 1963, RS codes employed a variable generator polynomial for encoding. This approach [peterson1972error](@cite) differed from the prevalent BCH scheme (used here), which utilizes a fixed generator polynomial. Consequently, these original RS codes weren't strictly categorized as BCH codes. Furthermore, depending on the chosen evaluation points, they might not even qualify as cyclic codes.
"""

function generator_polynomial(rs::ReedSolomon)
    GF2ʳ, a = finite_field(2, rs.m, "a")
    P, x = GF2ʳ[:x]
    gx = x - a ^ 1
    for i in 2:2 * rs.t
        gx *= (x - a ^ i)
    end
    return gx
end
