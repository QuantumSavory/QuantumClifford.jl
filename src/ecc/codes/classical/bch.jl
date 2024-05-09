"""The family of Bose–Chaudhuri–Hocquenghem (BCH) codes, as discovered in 1959 by Alexis Hocquenghem [hocquenghem1959codes](@cite), and independently in 1960 by Raj Chandra Bose and D.K. Ray-Chaudhuri [bose1960class](@cite).

The binary parity check matrix can be obtained from the following matrix over GF(2) field elements:

```
1		(α¹)¹			(α¹)²			(α¹)³			...		(α¹)ⁿ ⁻ ¹
1		(α³)¹			(α³)²			(α³)³			...		(α³)ⁿ ⁻ ¹
1		(α⁵)¹			(α⁵)²			(α⁵)³			...		(α⁵)ⁿ ⁻ ¹
.		  .			  .			  .			...		   .
.		  .			  .			  .			...		   .
.		  .			  .			  .			...		   .
1		(α²ᵗ ⁻ ¹)¹		(α²ᵗ ⁻ ¹)²		(α²ᵗ ⁻ ¹)³		...		(α²ᵗ ⁻ ¹)ⁿ ⁻ ¹
```

The entries of the matrix are in GF(2ᵐ). Each element in GF(2ᵐ) can be represented by an `m`-tuple (a binary column vector of length `m`). If each entry of `H` is replaced by its corresponding `m`-tuple, we obtain a binary parity check matrix for the code.

The BCH code is cyclic as its generator polynomial, `g(x)` divides `xⁿ - 1`, so `mod (xⁿ - 1, g(x)) = 0`.

You might be interested in consulting [bose1960further](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/q-ary_bch).
"""

abstract type AbstractPolynomialCode <: ClassicalCode end

"""
`BCH(m, t)`
- `m`: The positive integer defining the degree of the finite (Galois) field, GF(2ᵐ).
- `t`: The positive integer specifying the number of correctable errors.
"""
struct BCH <: AbstractPolynomialCode
    m::Int 
    t::Int 
    function BCH(m, t)
        if m < 3 || t < 0 || t >= 2 ^ (m - 1)
            throw(ArgumentError("Invalid parameters: `m` and `t` must be positive. Additionally, ensure `m ≥ 3` and `t < 2ᵐ ⁻ ¹` to obtain a valid code."))
        end
        new(m, t)
    end
end

"""
Generator Polynomial of BCH Codes

This function calculates the generator polynomial `g(x)` of a `t`-bit error-correcting BCH code of binary length `n = 2ᵐ - 1`. The binary code is derived from a code over the finite Galois field GF(2).

`generator_polynomial(BCH(m, t))`

- `m`: The positive integer defining the degree of the finite (Galois) field, GF(2ᵐ).
- `t`: The positive integer specifying the number of correctable errors.

Description:

The generator polynomial `g(x)` is the fundamental polynomial used for encoding and decoding BCH codes. It has the following properties:

1. Roots: It has `α`, `α²`, `α³`, ..., `α²ᵗ` as its roots, where `α` is a primitive element of the Galois Field GF(2ᵐ).
2. Error Correction: A BCH code with generator polynomial `g(x)` can correct up to `t` errors in a codeword of length `2ᵐ - 1`.
3. Minimal Polynomials: `g(x)` is the least common multiple (LCM) of the minimal polynomials `φᵢ(x)` of `αⁱ` for `i = 1` to `2ᵗ`.

Minimal Polynomial:

- The minimal polynomial of a field element `α` in GF(2ᵐ) is the polynomial of the lowest degree over GF(2) that has `α` as a root. It represents the simplest polynomial relationship between `α` and the elements of GF(2).

Least Common Multiple (LCM):

- The LCM of two or more polynomials `fᵢ(x)` is the polynomial with the lowest degree that is a multiple of all `fᵢ(x)`. It ensures that `g(x)` has all the roots of `φᵢ(x)` for `i = 1` to `2ᵗ`.
"""
function generator_polynomial(b::BCH)
    GF2ʳ, a = finite_field(2, b.m, "a")
    GF2x, x = GF(2)["x"]
    minimal_poly = FqPolyRingElem[]
    for i in 1:(2 * b.t - 1)
        if i % 2 != 0
            push!(minimal_poly, minpoly(GF2x, a ^ i))
        end 
    end
    gx = lcm(minimal_poly)
    return gx
end

function parity_checks(b::BCH)
    GF2ʳ, a = finite_field(2, b.m, "a")
    HField = Matrix{FqFieldElem}(undef, b.t, 2 ^ b.m - 1)
    for i in 1:b.t
        for j in 1:2 ^ b.m - 1
            base = 2 * i - 1  
            HField[i, j] = (a ^ base) ^ (j - 1)
        end
    end
    H = Matrix{Bool}(undef, b.m * b.t, 2 ^ b.m - 1)
    for i in 1:b.t
        row_start = (i - 1) * b.m + 1
        row_end = row_start + b.m - 1
        for j in 1:2 ^ b.m - 1
            t_tuple = Bool[]
            for k in 0:b.m - 1
                push!(t_tuple, !is_zero(coeff(HField[i, j], k)))
            end 
            H[row_start:row_end, j] .=  vec(t_tuple')
        end
    end 
    return H
end

code_n(b::BCH) = 2 ^ b.m - 1
code_k(b::BCH) = 2 ^ b.m - 1 - degree(generator_polynomial(BCH(b.m, b.t)))
