"""The family of `[2ᵐ + 1 - 3, k, 2ᵐ ⁺ ¹ - 3 - k]` Shortened Maximum Distance Separable `(ShortenedMDS)` codes are constructed from the `[2ᵐ + 1, k, 2ᵐ ⁺ ¹ - k]` Extended, Augmented Reed-Solomon codes from the corresponding first `x - 1` columns of latter's parity-check matrix, using `j = 0`, and setting `α₀, α₁, α₂, ..., αₓ ₋ ₁` to  `α⁰, α¹, α², ..., αˣ ⁻ ¹` in the parity-check matrix.

Denoting the `x` field elements as `0, α⁰, α¹, α²,... αˣ ⁻ ¹`, the shortened field parity-check matrix (`HF`) is given as follows:

```
(α⁰)ʲ			(α¹)ʲ					(α²)ʲ				...		(αˣ ⁻ ¹)ʲ
(α⁰)ʲ ⁺ ¹		(α¹)ʲ ⁺ ¹				(α²)ʲ ⁺ ¹			...		(αˣ ⁻ ¹)ʲ ⁺ ¹
(α⁰)ʲ ⁺ ²		(α¹)ʲ ⁺ ²				(α²)ʲ ⁺ ²			...		(αˣ ⁻ ¹)ʲ ⁺ ²
(α⁰)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹	(α¹)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹			(α²)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		...		(αˣ ⁻ ¹)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹
	.			.					.			...			.
	.			.					.			...			.
	.			.					.			...			.
(α⁰)ʲ ⁺ ˣ ⁻ ᵏ		(α¹)ʲ ⁺ ˣ ⁻ ᵏ	(α²)ʲ ⁺ ˣ ⁻ ᵏ			.			...		(αˣ ⁻ ¹)ʲ ⁺ ˣ ⁻ ᵏ
```

For significant coding gain, code length is typically restricted to less than 200 bits. Modified Dorsch decoder [dorsch1974decoding](@cite) is recommended for near maximum likelihood decoding.

Shortened MDS (`HF`) Matrix element expansion: 
    1. Row expansion: Each row of in the field parity-check matrix is replaced with an `m`-by-`m` field matrix defined over the base field `GF(2ᵐ)`.
    2. Column expansion: Consequently, the elements in each column of expanded field parity-check matrix are converted to binary representations by substituting powers of a primitive element (`α`) in the Galois Field `GF(2ᵐ)` with their corresponding `m`-tuples over the Boolean Field `GF(2)`.

You might be interested in consulting [tomlinson2017error](@cite), [macwilliams1977theory](@cite), [peterson1972error](@cite), [seroussi1986mds](@cite), [dur1987automorphism](@cite) and [Extending MDS Codes](https://www.unb.ca/faculty-staff/directory/_resources/pdf/sase/alderson/mds-codes.pdf) as well.
"""
abstract type AbstractPolynomialCode <: ClassicalCode end

struct ShortenedMDS <: AbstractPolynomialCode
    m::Int
    t::Int

    function ShortenedMDS(m, t)
        if m < 3 || t < 0 || t >= 2 ^ (m - 1) 
            throw(ArgumentError("Invalid parameters: m and t must be non-negative. Also, m > 3 and t < 2 ^ (m - 1) in order to obtain a valid code."))
        end
        new(m, t)
    end
end

"""
`parity_checks(ShortenedMDS(m, t))`
- `m`: The positive integer defining the degree of the finite (Galois) field, `GF(2ᵐ)`.
- `t`: The positive integer specifying the number of correctable errors.
"""
function parity_checks(rs::ShortenedMDS)
    GF2ʳ, a = finite_field(2, rs.m, "a")
    s_symbols = 3 # 3-level quantization. 
    x = 2 ^ rs.m + 1 - s_symbols
    k = 2 ^ rs.m - 1 - 2 * rs.t
    HField = Matrix{FqFieldElem}(undef, x - k + 1, x)
    for j in 1:x
        HField[1, j] = a ^ 0
    end
    for i in 1: x - k + 1
        HField[i, 1] = a ^ 0
    end
    for i in 2:x - k + 1
        for j in 2:x
            HField[i, j] = (a ^ (j - 1)) ^ (i - 2)
        end
    end
    HSeed = vcat(HField[1:1, :], HField[3:end, :])
    HTemp2 = Matrix{FqFieldElem}(undef, rs.m, x)
    HFieldExpanded = Matrix{FqFieldElem}(undef, rs.m * k, x)
    g = 1
    while g <= rs.m * k
        for i in 1:x - k
            for p in 1:rs.m
                HTemp2[p:p, :] = reshape(HSeed[i, :].*a ^ (p - 1) , 1, :)
            end
        if g > rs.m * k
           break
        end
        HFieldExpanded[g:g + rs.m - 1, :] .=  HTemp2
        g = g + rs.m
        end
    end
    H = Matrix{Bool}(undef, rs.m * k, rs.m * x)
    for i in 1:rs.m * k
        for j in 1:x
            col_start = (j - 1) * rs.m + 1
            col_end = col_start + rs.m - 1
            t_tuple = Bool[]
            for k in 0:rs.m - 1
                push!(t_tuple, !is_zero(coeff(HFieldExpanded[i, j], k)))
            end 
            H[i, col_start:col_end] .=  vec(t_tuple)
        end
    end
    return H
end

code_n(rs::ShortenedMDS) = (2 ^ rs.m + 1 - 3) * rs.m

code_k(rs::ShortenedMDS) = (2 ^ rs.m - 1 - 2 * rs.t) * rs.m
