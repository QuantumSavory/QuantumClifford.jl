"""The family of Bose–Chaudhuri–Hocquenghem (BCH) codes, as discovered in 1959 by Alexis Hocquenghem [hocquenghem1959codes](@cite), and independently in 1960 by Raj Chandra Bose and D.K. Ray-Chaudhuri [bose1960class](@cite).

The binary parity check matrix can be obtained from the following reduced matrix over `GF(2)` field elements:

```
1		(α¹)¹			(α¹)²			(α¹)³			...		(α¹)ⁿ ⁻ ¹
1		(α³)¹			(α³)²			(α³)³			...		(α³)ⁿ ⁻ ¹
1		(α⁵)¹			(α⁵)²			(α⁵)³			...		(α⁵)ⁿ ⁻ ¹
.		  .			  .			  .			...		   .
.		  .			  .			  .			...		   .
.		  .			  .			  .			...		   .
1		(α²ᵗ ⁻ ¹)¹		(α²ᵗ ⁻ ¹)²		(α²ᵗ ⁻ ¹)³		...		(α²ᵗ ⁻ ¹)ⁿ ⁻ ¹
```

The entries of the matrix are in `GF(qᵐ)`. Each element in `GF(qᵐ)` can be represented by an `m`-tuple (a binary column vector of length `m`). If each entry of `H` is replaced by its corresponding `m`-tuple, we obtain a binary parity check matrix for the code.

BCH codes come in two main categories depending on the value of `q` in the finite Galois field `GF(qᵐ)`:

1. Binary BCH Codes `(q = 2)`: Operate over a binary field (0, 1), making them efficient for binary communication.
2. Non-binary BCH Codes `(q ≠ 2)`: Utilize a finite Galois field where `q` is any prime number besides 2. Example: Reed-Solomon codes.

The BCH code is cyclic as its generator polynomial, `g(x)` divides `xⁿ - 1`, so `mod (xⁿ - 1, g(x)) = 0`.

You might be interested in consulting [bose1960further](@cite) and [error2024lin](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/q-ary_bch).
"""

abstract type AbstractPolynomialCode <: ClassicalCode end

"""
`BCH(q, m, t)`
- `q`: The positive prime number that determines the base size of the finite (Galois) field, GF(qᵐ). 
- `m`: The positive integer defining the degree of the finite (Galois) field, GF(qᵐ).
- `t`: The positive integer specifying the number of correctable errors.
"""
struct BCH <: AbstractPolynomialCode
    q::Int 
    m::Int 
    t::Int
    function BCH(q, m, t)
        if q < 2 || !isprime(q) || m < 3 || t < 0 || t >= q ^ (m - 1)
            throw(ArgumentError("Invalid parameters: `q`, `m` and `t` must be positive. Additionally, ensure that `q` is a prime number, `m ≥ 3` and `t < qᵐ ⁻ ¹` to obtain a valid code."))
        end
        new(q, m, t)
    end
end

"""
Generator Polynomial of BCH Codes

This function calculates the generator polynomial `g(x)` of a `t`-bit error-correcting BCH code of binary length `n = qᵐ - 1`. The binary code is derived from a code over the finite Galois field `GF(2)`.

`generator_polynomial(BCH(q, m, t))`
- `q`: The positive prime number that determines the base size of the finite (Galois) field, GF(qᵐ). 
- `m`: The positive integer defining the degree of the finite (Galois) field, GF(qᵐ).
- `t`: The positive integer specifying the number of correctable errors.

The generator polynomial `g(x)` is the fundamental polynomial used for encoding and decoding BCH codes. It has the following properties:

1. Roots: It has `α`, `α²`, `α³`, ..., `α²ᵗ` as its roots, where `α` is a primitive element of the Galois Field `GF(qᵐ)`.
2. Error Correction: A BCH code with generator polynomial `g(x)` can correct up to `t` errors in a codeword of length `qᵐ - 1`.
3. Minimal Polynomials: `g(x)` is the least common multiple (LCM) of the minimal polynomials `φᵢ(x)` of `αⁱ` for `i` from `1` to `2ᵗ`.

Useful definitions and background:

Minimal Polynomial: The minimal polynomial of a field element `α` in `GF(qᵐ)` is the polynomial of the lowest degree over `GF(2)` that has `α` as a root.

Least Common Multiple (LCM): The LCM of two or more polynomials `fᵢ(x)` is the polynomial with the lowest degree that is a multiple of all `fᵢ(x)`. It ensures that `g(x)` has all the roots of `φᵢ(x)` for `i = 1` to `2ᵗ`.

Conway polynomial: The finite Galois field `GF(qᵐ)` can have multiple distinct primitive polynomials of the same degree due to existence of several irreducible polynomials of that degree, each generating the field through different roots. Nemo.jl uses [Conway polynomial](https://en.wikipedia.org/wiki/Conway_polynomial_(finite_fields)), a standard way to represent the primitive polynomial for finite Galois fields `GF(pᵐ)` of degree `m`, where `p` is a prime number.
"""
function generator_polynomial(b::BCH)
    GFqʳ, a = finite_field(b.q, b.m, "a")
    GFqx, x = GF(b.q)["x"]
    minimal_poly = FqPolyRingElem[]
    for i in 1:2 * b.t
        if i % 2 != 0
            push!(minimal_poly, minpoly(GFqx, a ^ i))
        end 
    end
    gx = lcm(minimal_poly)
    return gx
end

function parity_checks(b::BCH)
    GFqʳ, a = finite_field(b.q, b.m, "a")
    n = b.q ^ b.m - 1
    HField = Matrix{FqFieldElem}(undef, b.t, n)
    for i in 1:b.t
        for j in 1:n
            base = 2 * i - 1  
            HField[i, j] = (a ^ base) ^ (j - 1)
        end
    end
    H = Matrix{Bool}(undef, b.m * b.t, n)
    for i in 1:b.t
        row_start = (i - 1) * b.m + 1
        row_end = row_start + b.m - 1
        for j in 1:n
            t_tuple = Bool[]
            for k in 0:b.m - 1
                push!(t_tuple, !is_zero(coeff(HField[i, j], k)))
            end 
            H[row_start:row_end, j] .=  vec(t_tuple')
        end
    end 
    return H
end

code_n(b::BCH) = b.q ^ b.m - 1
code_k(b::BCH) = b.q ^ b.m - 1 - degree(generator_polynomial(BCH(b.q, b.m, b.t)))
