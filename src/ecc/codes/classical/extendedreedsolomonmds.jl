"""Extended Reed-Solomon Maximum Distance Separable `(ExtendedReedSolomonMDS)` codes are constructed from the Galois Field `GF(2ᵐ)`. These codes are extensions of the family of Reed-Solomon codes, as discovered by Reed and Solomon in their 1960 paper [reed1960polynomial](@cite). 

These codes possess a code length `(n)` of `2ᵐ + 1` and are Maximum Distance Separable `(MDS)` codes with parameters `[[2ᵐ + 1, k, 2ᵐ⁺¹ − k]]`. These ExtendedReedSolomonMDS codes are not binary codes but frequently are used with `x = 2ᵐ`, and so there is a mapping of residue classes of a primitive polynomial with binary coefficients and each element of `GF(2ᵐ)` is represented as a binary `m`-tuple. Denoting the `x` field elements as `0, α⁰, α¹, α²,... αˣ ⁻ ¹`, the shortened field parity-check matrix (`HF`) is given as follows:

```
(α⁰)ʲ			(α¹)ʲ			(α²)ʲ				...		(αˣ ⁻ ¹)ʲ
(α⁰)ʲ ⁺ ¹		(α¹)ʲ ⁺ ¹		(α²)ʲ ⁺ ¹			...		(αˣ ⁻ ¹)ʲ ⁺ ¹
(α⁰)ʲ ⁺ ²		(α¹)ʲ ⁺ ²		(α²)ʲ ⁺ ²			...		(αˣ ⁻ ¹)ʲ ⁺ ²
(α⁰)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹	(α¹)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹	(α²)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		...		(αˣ ⁻ ¹)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹
	.			.			.			...			.
	.			.			.			...			.
	.			.			.			...			.
(α⁰)ʲ ⁺ ˣ ⁻ ᵏ		(α¹)ʲ ⁺ ˣ ⁻ ᵏ	(α²)ʲ ⁺ ˣ ⁻ ᵏ		...		(αˣ ⁻ ¹)ʲ ⁺ ˣ ⁻ ᵏ
```

You might be interested in consulting [tomlinson2017error](@cite), [macwilliams1977theory](@cite), [peterson1972error](@cite), [seroussi1986mds](@cite), [dur1987automorphism](@cite) and [Extending MDS Codes](https://www.unb.ca/faculty-staff/directory/_resources/pdf/sase/alderson/mds-codes.pdf) as well.

Extended RS codes may be constructed using any Galois Field `GF(x)`, resulting in parameters [[x + 1, k, x + 2 − k]] [macwilliams1977theory](@cite). The ECC Zoo has an [entry for Extended Generalized Reed-Solomon codes](https://errorcorrectionzoo.org/c/extended_reed_solomon).
"""

abstract type AbstractPolynomialCode <: ClassicalCode end

struct ExtendedReedSolomonMDS <: AbstractPolynomialCode
    m::Int
    t::Int

    function ExtendedReedSolomonMDS(m, t)
        if m < 3 || t < 0 || t >= 2 ^ (m - 1) 
            throw(ArgumentError("Invalid parameters: m and t must be non-negative. Also, m > 3 and t < 2 ^ (m - 1) in order to obtain a valid code."))
        end
        new(m, t)
    end
end

"""
`parity_checks(ExtendedReedSolomonMDS(m, t))`
- `m`: The positive integer defining the degree of the finite (Galois) field, `GF(2ᵐ)`.
- `t`: The positive integer specifying the number of correctable errors.

This function applies Extended Reed-Solomon Maximum Distance Separable `(ExtendedReedSolomonMDS)` codes for binary transmission using soft decisions (see section 7.3)[tomlinson2017error](@cite). For significant coding gain, code length is typically restricted to less than 200 bits. Modified Dorsch decoder [dorsch1974decoding](@cite) is recommended for near maximum likelihood decoding.

Challenges of Standard RS Codes: While efficient as MDS codes, standard RS codes are not ideal for binary channels. As demonstrated in the results (see section 7.2)[tomlinson2017error](@cite), their performance suffers due to a mismatch between the code structure (symbol-based) and the channel (binary). A single bit error can lead to a symbol error, negating the code's benefits.

Improved Binary Codes through Concatenation: This method enhances RS codes for binary channels through code concatenation. It adds a single overall binary parity check to each `m`-tuple representing a symbol. This approach transforms the original RS code `[[n, k, n - k - 1]]` into a new binary code with parameters `[[n[m + 1], k * m, 2[n - k -1]]]`. The resulting binary code boasts a minimum symbol weight of 2, effectively doubling the minimum Hamming distance compared to the original RS code.

Augmented Extended RS Codes: Constructed from Galois Field `GF(2ᵐ)`. Length: `2ᵐ + 1` (Maximum Distance Separable (MDS) codes). Parameters: `[[2ᵐ + 1, k, 2ᵐ ⁺ ¹ - k]]`. Generalization: Applicable to any Galois Field `GF(x)` with parameters `[[x + 1, k, x + 2 - k]]`.

Field Parity-Check Matrix Properties:

```
(α₀)ʲ				(α₁)ʲ				(α₂)ʲ			...		(αₓ₋₂)ʲ				1	0
(α₀)ʲ ⁺ ¹			(α₁)ʲ ⁺ ¹			(α₂)ʲ ⁺ ¹		...		(αₓ₋₂)ʲ ⁺ ¹			0	0
(α₀)ʲ ⁺ ²			(α₁)ʲ ⁺ ²			(α₂)ʲ ⁺ ²		...		(αₓ₋₂)ʲ ⁺ ²			0	0
(α₀)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		(α₁)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		(α₂)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹	...		(αₓ₋₂)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		0	0
	.				.				.		...			.			.	.
	.				.				.		...			.			.	.
	.				.				.		...			.			.	.
(α₀)ʲ ⁺ ˣ ⁻ ᵏ			(α₁)ʲ ⁺ ˣ ⁻ ᵏ			(α₂)ʲ ⁺ ˣ ⁻ ᵏ		...		(αₓ₋₂)ʲ ⁺ ˣ ⁻ ᵏ		0	1
```

The matrix has `x - k + 1` rows corresponding to the code's parity symbols. Any `x - k + 1` columns form a Vandermonde matrix (non-singular). This ensures correction of up to `x - k + 1` symbol erasures in a codeword. We can re-arrange the columns of this matrix in any desired order. Any set of `s` symbols within a codeword can be designated as parity symbols and permanently removed. This important property leads to construction of Shortened MDS codes.

Shortened MDS Codes: Corresponding columns of the field parity-check matrix can be deleted to form a shortened `[[2ᵐ + 1 - s, k, 2ᵐ ⁺ ¹ - s - k]]` MDS code. This is an important property of MDS codes, particularly for their practical realisation in the form of augmented, extended RS codes because it enables efficient implementation in applications such as incremental redundancy systems, and network coding. The 3-level quantization of the received channel bits is utilized meaning 3 symbols are deleted. The Fig. 7.2 [tomlinson2017error](@cite) shows that with 3-level quantization, there is an improvement over the binary-transmission with hard decisions for Reed-Solomon coding. The designed distance for binary expanded parity check matrix remains same as symbol based parity check matrix. According to [macwilliams1977theory](@cite), changing the basis `j` can increase the designed distance `(dmin)` of the resulting binary code.

Cyclic Code Construction: Using the first `x - 1` columns of the field parity-check matrix, using `j = 0`, and setting `α₀, α₁, α₂, ..., αₓ ₋ ₁` to  `α⁰, α¹, α², ..., αˣ ⁻ ¹` in the parity-check matrix are set equal to the powers of a primitive element α of the Galois Field `GF(x)`, a cyclic code can be constructed for efficient encoding and decoding.

Shortened MDS (`HF`) Matrix element expansion: 
    1. Row expansion: Each row of in the field parity-check matrix is replaced with an `m`-by-`m` field matrix defined over the base field `GF(2ᵐ)`.
    2. Column expansion: Consequently, the elements in each column of expanded field parity-check matrix are converted to binary representations by substituting powers of a primitive element (`α`) in the Galois Field `GF(2ᵐ)` with their corresponding `m`-tuples over the Boolean Field `GF(2)`.
"""
function parity_checks(rs::ExtendedReedSolomonMDS)
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

code_n(rs::ExtendedReedSolomonMDS) = (2 ^ rs.m + 1 - 3) * rs.m
code_k(rs::ExtendedReedSolomonMDS) = (2 ^ rs.m - 1 - 2 * rs.t) * rs.m
