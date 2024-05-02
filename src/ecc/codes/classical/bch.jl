"""The family of Bose–Chaudhuri–Hocquenghem (BCH) codes, as discovered in 1959 by Alexis Hocquenghem [hocquenghem1959codes](@cite), and independently in 1960 by Raj Chandra Bose and D.K. Ray-Chaudhuri [bose1960class](@cite).

The binary parity check matrix can be obtained from the following matrix over field elements, after each field element is expressed as a binary column vector over GF(2).

```
1            α¹                       α²                     α³                   ...              αⁿ ⁻ ¹
1            (α ^ 3)¹                 (α ^ 3)²               (α ^ 3)³             ...              (α ^ 3)ⁿ ⁻ ¹
1            (α ^ 5)¹                 (α ^ 5)²               (α ^ 5)³             ...              (α ^ 5)ⁿ ⁻ ¹
.                  .                         .                      .                                      .
.                  .                         .                      .                                      .
.                  .                         .                      .                                      .
1            (α ^ (2*t - 1))¹         (α ^ (2*t - 1))²       (α ^ (2*t - 1))³     ...              (α ^ (2*t - 1))ⁿ ⁻ ¹

```

BCH code is cyclic code as its generator polynomial, `gx` divides `xⁿ - 1`, so `mod (xⁿ - 1, gx)` = 0 

Note: The entries of matrix over field elements are in GF(2ᵐ). Each element in GF(2ᵐ) can be represented by a m-tuple/binary column of length m over GF(2). If each entry of H is replaced by its corresponding m-tuple/binary column of length m over GF(2) arranged in column form, we obtain a binary parity check matrix for the code.

You might be interested in consulting [bose1960further](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/q-ary_bch).
"""

abstract type AbstractPolynomialCode <: ClassicalCode end

"""
BCH(m, t):
- `m` (Integer): The positive integer defining the degree of the finite (Galois) field, GF(2ᵐ).
- `t` (Integer): The positive integer specifying the number of correctable errors (`t`).
"""
struct BCH <: AbstractPolynomialCode
    m::Int 
    t::Int 
    function BCH(m, t)
        if m < 3 || m > 10 || t < 0 || t >= 2^(m - 1)
            throw(ArgumentError("Invalid parameters: 'm' and 't' must be positive. Additionally, 't' < 2ᵐ ⁻ ¹ to obtain a valid code and to tractable."))
        end
        new(m, t)
    end
end

"""
Generator Polynomial of BCH Codes

This function calculates the generator polynomial `g(x)` of a t-bit error-correcting BCH code of length `(n)` of `2ᵐ - 1` over the finite Galois field GF(2).

generator_polynomia(BCH(m, t)):
- `m` (Integer): The positive integer defining the degree of the finite (Galois) field, GF(2ᵐ).
- `t` (Integer): The positive integer specifying the number of correctable errors (`t`).

Description:

The generator polynomial `g(x)` is the fundamental polynomial used for encoding and decoding BCH codes. It has the following properties:

1. Roots: It has `α`, `α²`, `α³`, ..., `α²ᵗ` as its roots, where `α` is a primitive element of the Galois Field GF(2ᵐ).
2. Error Correction: A BCH code with generator polynomial `g(x)` can correct up to `t` errors in a codeword of length `2ᵐ - 1`.
3. Minimal Polynomials: `g(x)` is the least common multiple (LCM) of the minimal polynomials `φ_i(x)` of `αⁱ` for `i = 1` to `2ᵗ`.

Minimal Polynomial:

- The minimal polynomial of a field element `α` in GF(2ᵐ) is the polynomial of the lowest degree over GF(2) that has `α` as a root. It represents the simplest polynomial relationship between `α` and the elements of GF(2).

Least Common Multiple (LCM):

- The LCM of two or more polynomials `f_i(x)` is the polynomial with the lowest degree that is a multiple of all `f_i(x)`. It ensures that `g(x)` has all the roots of `φ_i(x)` for `i = 1` to `2ᵗ`.
"""
function generator_polynomial(rs::BCH)
    GF2ͬ, a = finite_field(2, rs.m, "a")
    GF2x, x = GF(2)["x"]
    minimal_poly = FqPolyRingElem[]
    for i in 1:(2*rs.t - 1)
        if i % 2 != 0
            push!(minimal_poly, minpoly(GF2x, a^i))
        end 
    end
    gx = lcm(minimal_poly)
    return gx
end

function parity_checks(rs::BCH)
    GF2ͬ, a = finite_field(2, rs.m, "a")
    HField = Matrix{FqFieldElem}(undef, rs.t, 2^rs.m - 1)
    for i in 1:rs.t
        for j in 1:2^rs.m - 1
            base = 2 * i - 1  
            HField[i, j] = (a ^ base) ^ (j - 1)
        end
    end
    H = Matrix{Bool}(undef, rs.m*rs.t, 2^rs.m - 1)
    for i in 1:rs.t
        row_start = (i - 1) * rs.m + 1
        row_end = row_start + rs.m - 1
        for j in 1:2^rs.m - 1
            t_tuple = Bool[]
            for k in 0:rs.m - 1
                push!(t_tuple, !is_zero(coeff(HField[i, j], k)))
            end 
            H[row_start:row_end, j] .=  vec(t_tuple')
        end
    end 
    return H
end