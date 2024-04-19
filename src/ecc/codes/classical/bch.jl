"""The family of Bose–Chaudhuri–Hocquenghem (BCH) codes, as discovered in 1959 by Alexis Hocquenghem [hocquenghem1959codes](@cite), and independently in 1960 by Raj Chandra Bose and D.K. Ray-Chaudhuri [bose1960class](@cite).

The binary parity check matrix can be obtained from the following matrix over field elements, after each field element is expressed as a binary column vector over GF(2).

1            a^1             a^2                a^3               ...          a^(n - 1)
1         (a^3)^1          (a^3)^2            (a^3)^3             ...       (a^3)^(n - 1)
1         (a^5)^1          (a^5)^2            (a^5)^3             ...       (a^5)^(n - 1)
.            .               .                   .                                  .
.            .               .                   .                                  .
.            .               .                   .                                  .
1       a^(2*t - 1)     a^(2*t - 1)^2      a^(2*t - 1)^3          ...      a^(2*t - 1)^(n -1) 

Note: The entries of matrix over field elements are in GF(2^m). Each element in GF(2^m) can be represented by a m-tuple/binary column of length m over GF(2). If each entry of H is replaced by its corresponding m-tuple/binary column of length m over GF(2) arranged in column form, we obtain a binary parity check matrix for the code.

You might be interested in consulting [bose1960further](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/q-ary_bch)
"""
struct BCH <: ClassicalCode
    n::Int 
    t::Int #Error correction capability; t bits can be corrected

    function BCH(n, t)
        if n < 6 || n > 500 || t < 0 || t > 2^(ceil(Int, log2(n + 1)) - 1)
            throw(ArgumentError("Invalid parameters: 'n' and 't' must be positive, and 'r' must be >= to 3. Additionally, 'n' is >= to 7 since n = 2ͬ - 1 and 't' < 2^(r - 1), to obtain a valid code and to tractable."))
        end
        new(n, t)
    end
end

"""
Generator Polynomial of BCH Codes

This function calculates the generator polynomial `g(x)` of a t-bit error-correcting BCH code of length `2^m - 1` over the finite Galois field GF(2).

Input Arguments:
- `m` (Integer): The positive integer defining the code length (`2^m - 1`).
- `t` (Integer): The positive integer specifying the number of correctable errors (`t`).

Description:

The generator polynomial `g(x)` is the fundamental polynomial used for encoding and decoding BCH codes. It has the following properties:

1. Roots: It has `α`, `α^2`, `α^3`, ..., `α^(2^t)` as its roots, where `α` is a primitive element of the Galois Field GF(2^m).
2. Error Correction: A BCH code with generator polynomial `g(x)` can correct up to `t` errors in a codeword of length `2^m - 1`.
3. Minimal Polynomials: `g(x)` is the least common multiple (LCM) of the minimal polynomials `φ_i(x)` of `α^i` for `i = 1` to `2^t`.

Minimal Polynomial:

- The minimal polynomial of a field element `α` in GF(2^m) is the polynomial of the lowest degree over GF(2) that has `α` as a root. It represents the simplest polynomial relationship between `α` and the elements of GF(2).

Least Common Multiple (LCM):

- The LCM of two or more polynomials `f_i(x)` is the polynomial with the lowest degree that is a multiple of all `f_i(x)`. It ensures that `g(x)` has all the roots of `φ_i(x)` for `i = 1` to `2^t`.
"""
function generator_polynomial(rs::BCH)
    r = ceil(Int, log2(rs.n + 1))
    GF2ͬ, a = finite_field(2, r, "a")
    GF2x, x = GF(2)["x"]
    minimal_poly = FqPolyRingElem[]
    for i in 1:(2*rs.t - 1)
        if i % 2 != 0
             minimal_poly = [minimal_poly; minpoly(GF2x, a^i)]
        end 
    end
    gx = lcm(minimal_poly)
    return gx
end

function parity_checks(rs::BCH)
    r = ceil(Int, log2(rs.n + 1))
    GF2ͬ, a = finite_field(2, r, "a")
    HField = Matrix{FqFieldElem}(undef, rs.t, rs.n)
    for i in 1:rs.t
        for j in 1:rs.n
            base = 2*i - 1  
            HField[i, j] = (a^base)^(j - 1)
        end
    end
    H = Matrix{Bool}(undef, r*rs.t, rs.n)
    for i in 1:rs.t
        row_start = (i - 1) * r + 1
        row_end = row_start + r - 1
        for j in 1:rs.n
            t_tuple = Bool[]
            for k in 0:r - 1
                t_tuple = [t_tuple; !is_zero(coeff(HField[i, j], k))]
            end 
            H[row_start:row_end, j] .=  vec(t_tuple')
        end
    end 
    return H
end
