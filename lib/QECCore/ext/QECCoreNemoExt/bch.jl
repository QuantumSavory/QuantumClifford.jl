"""
Generator Polynomial of BCH Codes

This function calculates the generator polynomial `g(x)` of a `t`-bit error-correcting
BCH code of binary length `n = 2ᵐ - 1`. The binary code is derived from a code over the
finite Galois field `GF(2ᵐ)`.

The generator polynomial `g(x)` is the fundamental polynomial used for encoding and
decoding BCH codes. It has the following properties:

- Roots: It has `α`, `α²`, `α³`, ..., `α²ᵗ` as its roots, where `α` is a primitive element
of the Galois Field `GF(2ᵐ)`.
- Error Correction: A BCH code with generator polynomial `g(x)` can correct up to `t` errors
in a codeword of length `2ᵐ - 1`.
- Minimal Polynomials: `g(x)` is the least common multiple (LCM) of the minimal polynomials
`φᵢ(x)` of `αⁱ` for `i` from `1` to `2ᵗ`.

Useful definitions and background:

Minimal Polynomial: The minimal polynomial of a field element `α` in GF(2ᵐ) is the polynomial
of the lowest degree over `GF(2ᵐ)` that has `α` as a root.

Least Common Multiple (LCM): The LCM of two or more polynomials `fᵢ(x)` is the polynomial
with the lowest degree that is a multiple of all `fᵢ(x)`. It ensures that `g(x)` has all
the roots of `φᵢ(x)` for `i = 1` to `2ᵗ`.

Conway polynomial: The finite Galois field `GF(2ᵐ)` can have multiple distinct primitive
polynomials of the same degree due to existence of several irreducible polynomials of that
degree, each generating the field through different roots.

Nemo.jl uses [Conway polynomial](https://en.wikipedia.org/wiki/Conway_polynomial_(finite_fields)),
a standard way to represent the primitive polynomial for finite Galois fields `GF(pᵐ)` of degree
`m`, where `p` is a prime number.
"""
function generator_polynomial(b::BCH)
    GF2ʳ, a = finite_field(2, b.m, "a")
    GF2x, x = GF(2)["x"]
    minimal_poly = FqPolyRingElem[]
    for i in 1:2 * b.t
        if i % 2 != 0
            push!(minimal_poly, minpoly(GF2x, a ^ i))
        end
    end
    gx = lcm(minimal_poly)
    return gx
end

function parity_matrix(b::BCH)
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

code_k(b::BCH) = code_n(b) - degree(generator_polynomial(b))
