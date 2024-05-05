"""The family of Reed-Solomon codes, as discovered by Reed and Solomon in their 1960 paper [reed1960polynomial](@cite). 

Reed Solomon codes are maximum distance separable (MDS) codes and have the highest possible minimum Hamming distance. The codes have symbols from F_q with parameters [[x - 1, k, x - k]].

They are not binary codes but frequently are used with x = 2ᵐ, and so there is a mapping of residue classes of a primitive polynomial with binary coefficients and each element of GF(2ᵐ) is represented as a binary m-tuple. Denoting the x field elements as 0, α⁰, α¹, α²,... αˣ ⁻ ¹, the shortened Field parity-check matrix (`HSeed`) is given by

```
(α⁰)ʲ			(α¹)ʲ			(α²)ʲ				...		(αˣ ⁻ ¹)ʲ
(α⁰)ʲ ⁺ ¹		(α¹)ʲ ⁺ ¹		(α²)ʲ ⁺ ¹			...		(αˣ ⁻ ¹)ʲ ⁺ ¹
(α⁰)ʲ ⁺ ²		(α¹)ʲ ⁺ ²		(α²)ʲ ⁺ ²			...		(αˣ ⁻ ¹)ʲ ⁺ ²
(α⁰)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹	(α¹)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹	(α²)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		...		(αˣ ⁻ ¹)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹
	.			.			.			...			.
	.			.			.			...			.
	.			.			.			...			.
(α⁰)ʲ ⁺ ˣ ⁻ ᵏ		(α¹)ʲ ⁺ ˣ ⁻ ᵏ		(α²)ʲ ⁺ ˣ ⁻ ᵏ			...		(αˣ ⁻ ¹)ʲ ⁺ ˣ ⁻ ᵏ
```

You might be interested in consulting [geisel1990tutorial](@cite), [wicker1999reed](@cite), [sklar2001reed](@cite), [berlekamp1978readable](@cite), [tomlinson2017error](@cite), [macwilliams1977theory](@cite), [https://youtu.be/K26Ssr8H3ec?si=QOeohq_6I0Oyd8qu] and [peterson1972error](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/reed_solomon)
"""

abstract type AbstractPolynomialCode <: ClassicalCode end

struct ReedSolomon <: AbstractPolynomialCode
    m::Int
    t::Int

    function ReedSolomon(m, t)
        if m < 3 || t < 0 || t >= 2 ^ (m - 1) 
            throw(ArgumentError("Invalid parameters: m and t must be non-negative. Also, m > 3 and t < 2 ^ (m - 1) in order to obtain a valid code and to remain tractable."))
        end
        new(m, t)
    end
end

"""
generator_polynomial(ReedSolomon(m, t)):
- `m` (Integer): The positive integer defining the degree of the finite (Galois) field, GF(2ᵐ).
- `t` (Integer): The positive integer specifying the number of correctable errors (`t`).

The generator polynomial for an RS code takes the following form:

```
g(X) = g₀ + g₁X¹ + g₂X² + ... + g₂ₜ₋₁X²ᵗ⁻¹ + X²ᵗ
```

where:
- `X` is the indeterminate variable.
- `gᵢ` are the coefficients of the polynomial.
- `t` is the number of correctable symbol errors.

We describe the generator polynomial in terms of its 2 * `t`  = `n` - `k` roots, as follows:

``` 
g(X) = (X - α¹)(X - α²)(X - α³) ... (X - α²ᵗ)
```

Degree and Parity Symbols:

- The degree of the generator polynomial is equal to 2 * `t`, which is also the number of parity symbols added to the original data (`k` symbols) to create a codeword of length `n` (`n` = `k` + 2 * `t`).

Roots of the Generator Polynomial:

- The generator polynomial has 2 * `t` distinct roots, designated as `α¹, α², ... , α²ᵗ`. These roots are chosen from a Finite Galois Field.
- Any power of α can be used as the starting root, not necessarily `α¹` itself.

Notes:

- Only in this construction scheme for defining the fixed generator polynomial, `g(x)`, RS codes are a subset of the Bose, Chaudhuri, and Hocquenghem (BCH) codes; hence, this relationship between the degree of the generator polynomial and the number of parity symbols holds, just as for BCH codes where degree of BCH generator polynomial, `degree(g(x)) == n - k`. This is proved via RS test suite as well.

- Prior to 1963, RS codes employed a variable generator polynomial for encoding. This approach [peterson1972error](@cite) differed from the prevalent BCH scheme (used here), which utilizes a fixed generator polynomial. Consequently, these original RS codes weren't strictly categorized as BCH codes. Furthermore, depending on the chosen evaluation points, they might not even qualify as cyclic codes.

"""
function generator_polynomial(rs::ReedSolomon)
    GF2ʳ, a = finite_field(2, rs.m, "a")
    P, x = GF2ʳ[:x]
    gx = x - a ^ 1
    for i in 2:2*rs.t
        gx *= (x - a ^ i)
    end
    return gx
end

"""
This function applies Reed-Solomon (RS) codes with soft decision decoding for binary transmission channels.

Challenges of Standard RS Codes:

- While efficient as MDS codes, standard RS codes are not ideal for binary channels.
- As demonstrated in the results (see section 7.2)[tomlinson2017error](@cite), their performance suffers due to a mismatch between the code structure (symbol-based) and the channel (binary).
- A single bit error can lead to a symbol error, negating the code's benefits.

Improved Binary Codes through Concatenation:

- This method enhances RS codes for binary channels through code concatenation.
- It adds a single overall binary parity check to each m-tuple representing a symbol.
- This approach transforms the original RS code [[n, k, n - k - 1]] into a new binary code with parameters [[n[m + 1], k * m, 2[n - k -1]]].
- The resulting binary code boasts a minimum symbol weight of 2, effectively doubling the minimum Hamming distance compared to the original RS code.

Key Points:

- Uses unquantized soft decision decoding for improved performance.
- Modified Dorsch decoder is recommended for near maximum likelihood decoding. 
- Code length limitations: For significant coding gain, code length is typically restricted to less than 200 bits.

Augmented Extended RS Codes:

- Constructed from Galois Field GF(2ᵐ).
- Length: 2ᵐ + 1 (Maximum Distance Separable (MDS) codes).
- Parameters: [[2ᵐ + 1, k, 2ᵐ ⁺ ¹ - k]].
- Generalization: Applicable to any Galois Field GF(q) with parameters [[x + 1, k, x + 2 - k]].

Field Parity-Check Matrix (`HField`) Properties:

```
(α₀)ʲ				(α₁)ʲ				(α₂)ʲ			...		(αₓ₋₂)ʲ				1	0
(α₀)ʲ ⁺ ¹			(α₁)ʲ ⁺ ¹			(α₂)ʲ ⁺ ¹		...		(αₓ₋₂)ʲ ⁺ ¹			0	0
(α₀)ʲ ⁺ ²			(α₁)ʲ ⁺ ²			(α₂)ʲ ⁺ ²		...		(αₓ₋₂)ʲ ⁺ ²			0	0
(α₀)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		(α₁)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		(α₂)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹	...		(αₓ₋₂)ʲ ⁺ ˣ ⁻ ᵏ ⁻ ¹		0	0
	.				.				.		...			.			.	.
	.				.				.		...			.			.	.
	.				.				.		...			.			.	.
(α₀)ʲ ⁺ ˣ ⁻ ᵏ			(α₁)ʲ ⁺ ˣ ⁻ ᵏ			(α₂)ʲ ⁺ ˣ ⁻ ᵏ		...		(αₓ₋₂)ʲ ⁺ ˣ ⁻ ᵏ			0	1
```

- The matrix has x - k + 1 rows corresponding to the code's parity symbols.
- Any x - k + 1 columns form a Vandermonde matrix (non-singular).
- This ensures correction of up to x - k + 1 symbol erasures in a codeword.
- Permutation: We can re-arrange the columns of the `HField` matrix in any desired order.
- Parity symbols (s) deletion: Any set of `s` symbols within a codeword can be designated as parity symbols and permanently removed. This important property leads to construction of Shortened MDS codes.

Shortened MDS Codes:

- Corresponding columns of the field parity-check matrix `HField` can be deleted to form a shortened [[2ᵐ + 1 - s, k, 2ᵐ ⁺ ¹ - s - k]] MDS code.
- This is an important property of MDS codes, particularly for their practical realisation in the form of augmented, extended RS codes because it enables efficient implementation in applications such as incremental redundancy systems, and network coding.
- 3 - level quantization of the received channel bits meaning 3 symbols deleted

Cyclic Code Construction:

- Using the first x - 1 columns of the field parity-check matrix (HField), using j = 0, and setting α₀, α₁, α₂, ..., αₓ ₋ ₁ to  α⁰, α¹, α², ..., αˣ ⁻ ¹ in the parity-check matrix are set equal to the powers of a primitive element α of the Galois Field GF(q), a cyclic code can be constructed for efficient encoding and decoding. The resulting matrix is represented by `HSeed`.

`HSeed` Matrix element expansion: 

    1. Row expansion: Each row of in the `HField` matrix is replaced with an `m`-by-`m` Field matrix defined over the base field GF (2ᵐ). This expansion is represented by `HFieldExpanded`.
    2. Column expansion: The elements in each column of `HFieldExpanded` matrix are converted to binary representations by substituting powers of a primitive element (`α`) in the Galois Field GF(2ᵐ) with their corresponding m-tuples over the Boolean/Binary Field GF(2).
"""
function parity_checks(rs::ReedSolomon)
    GF2ʳ, a = finite_field(2, rs.m, "a")
    s_symbols = 3 # 3-level quantization. 
    x = 2 ^ rs.m + 1 - s_symbols
    k = 2 ^ rs.m - 1 - 2 * rs.t
    HField = Matrix{FqFieldElem}(undef, x - k + 1, x)
    for j in 1:x
        HField[1, j] = a ^ 0
    end
    HTemp2 = Matrix{FqFieldElem}(undef, rs.m, x)
    for i in 1: x - k + 1
        HField[i, 1] = a ^ 0
    end
    for i in 2:x - k + 1
        for j in 2:x
            HField[i, j] = (a ^ (j - 1)) ^ (i - 2)
        end
    end
    HSeed = vcat(HField[1:1, :], HField[3:end, :])
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
