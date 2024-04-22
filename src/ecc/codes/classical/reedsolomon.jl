"""The family of Reed-Solomon codes, as discovered by Reed and Solomon in their 1960 paper [reed1960polynomial](@cite). 

Reed Solomon codes are maximum distance separable (MDS) codes and have the highest possible minimum Hamming distance. The codes have symbols from F_q with parameters [[q - 1, k, q - k]].

They are not binary codes but frequently are used with q = 2 ^ m, and so there is a mapping of residue classes of a primitive polynomial with binary coefficients and each element of GF(2 ^ m) is represented as a binary m-tuple. Denoting the q field elements as 0, α ^ 0, α ^ 1, α ^ 2,... α ^ q - 1,
the shortened Field parity-check matrix (`HSeed`) is given by

(α ^ 0) ^ (j)                  (α ^ 1) ^ (j)                    (α ^ 2) ^ (j)                    ...        (α ^ (q - 1)) ^ (j)            
(α ^ 0) ^ (j + 1)              (α ^ 1) ^ (j + 1)                (α ^ 2) ^ (j + 1)                ...        (α ^ (q - 1)) ^ (j + 1)            
(α ^ 0) ^ (j + 2)              (α ^ 1) ^ (j + 2)                (α ^ 2) ^ (j + 2)                ...        (α ^ (q - 1)) ^ (j + 2)           
        .                              .                                   .                     ...                   .                     
        .                              .                                   .                     ...                   .                    
        .                              .                                   .                     ...                   .                     
(α ^ 0) ^ (j + q - k - 1)      (α ^ 1) ^ (j + q - k - 1)        (α ^ 2) ^ (j + q - k - 1)        ...        (α ^ (q - 1)) ^ (j + q - k - 1)    
(α ^ 0) ^ (j + q - k)          (α ^ 1) ^ (j + q - k)            (α ^ 2) ^ (j + q - k)            ...        (α ^ (q - 1)) ^ (j + q - k)        

You might be interested in consulting [geisel1990tutorial](cite), [wicker1999reed](@cite), [sklar2001reed](@cite), [berlekamp1978readable](cite), [tomlinson2017error](@cite) and [https://youtu.be/K26Ssr8H3ec?si=QOeohq_6I0Oyd8qu](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/reed_solomon)
"""
struct ReedSolomon <: ClassicalCode
    n::Int
    k::Int

    function ReedSolomon(n, k)
        if n < 0 || k < 0 || n > 500
            throw(ArgumentError("Invalid parameters: n and k must be non-negative and n > 500 in order to obtain a valid code and to remain tractable"))
        end
        new(n, k)
    end
end

#section 3.3.2 [geisel1990tutorial](cite) for constructing custom message polynomial m(x)
function _message_polynomial_rs(k, a, an, x, positions)
    message = 0*x
    for pos in positions
        if pos <= k
            message += a^(an)*x^(pos)
        else
            throw(DomainError("Invalid bit positions: The number of bit positions [0th,..., kth] to assign coefficient values should be <= k"))
        end
    end
    return message
end

#section 3.3.2 [geisel1990tutorial](cite) for constructing parity_check polynomial ck(x) given m(x) and g(x)
function _parity_check_polynomial_rs(n, k, x, message, gx)
    x_nminusk = x^(n - k)
    ck = mod(x_nminusk*message, gx)
    return ck
end

#section 3.3.2 [geisel1990tutorial](cite) for constructing codeword polynomial c(x) given m(x) and ck(x)
function _codeword_polynomial_rs(message, ck)
    cx = message + ck
    return cx
end

function generator_polynomial(rs::ReedSolomon)
    r = ceil(Int, log2(rs.n + 1))
    t = div(rs.n - rs.k, 2)
    GF2ͬ, a = finite_field(2, r, "a")
    P, x = GF2ͬ[:x]
    pzeros = 2*t
    gx = x - a^pzeros
    for i in 1:(2*t - 1)
       gx *= (x - a^(pzeros + i))
    end
    return gx
end

# Reed-Solomon Codes and Binary Transmission with soft decisions
"""
This function applies Reed-Solomon (RS) codes with soft decision decoding for binary transmission channels.

Challenges of Standard RS Codes:

- While efficient as MDS codes, standard RS codes are not ideal for binary channels.
- As demonstrated in the results (see section 7.2)[tomlinson2017error](@cite), their performance suffers due to a mismatch between the code structure (symbol-based) and the channel (binary).
- A single bit error can lead to a symbol error, negating the code's benefits.

Improved Binary Codes through Concatenation:

- This method enhances RS codes for binary channels through code concatenation.
- It adds a single overall binary parity check to each m-tuple representing a symbol.
- This approach transforms the original RS code [[n, k, n - k - 1]] into a new binary code with parameters [[n[m + 1], k*m, 2[n - k -1]]].
- The resulting binary code boasts a minimum symbol weight of 2, effectively doubling the minimum Hamming distance compared to the original RS code.

Key Points:

- Uses unquantized soft decision decoding for improved performance.
- Modified Dorsch decoder is recommended for near maximum likelihood decoding. 
- Code length limitations: For significant coding gain, code length is typically restricted to less than 200 bits.

Augmented Extended RS Codes:

- Constructed from Galois Field GF(2 ^ m).
- Length: 2 ^ m + 1 (Maximum Distance Separable (MDS) codes).
- Parameters: [[2 ^(m) + 1, k, 2 ^ (m + 1) - k]].
- Generalization: Applicable to any Galois Field GF(q) with parameters [[q + 1, k, q + 2 - k]].

Field Parity-Check Matrix (`HField`) Properties:

(α_0) ^ (j)                  (α_1) ^ (j)                    (α_2) ^ (j)                   ...          (α_(q - 2)) ^ (j)                1   0 
(α_0) ^ (j + 1)              (α_1) ^ (j + 1)                (α_2) ^ (j + 1)               ...          (α_(q - 2)) ^ (j + 1)            0   0
(α_0) ^ (j + 2)              (α_1) ^ (j + 2)                (α_2) ^ (j + 2)               ...          (α_(q - 2)) ^ (j + 2)            0   0 
       .                            .                                .                    ...                   .                       .   .
       .                            .                                .                    ...                   .                       .   .
       .                            .                                .                    ...                   .                       .   .
(α_0) ^ (j + q - k - 1)      (α_1) ^ (j + q - k - 1)        (α_2) ^ (j + q - k - 1)       ...          (α_(q - 2)) ^ (j + q - k - 1)    0   0
(α_0) ^ (j + q - k)          (α_1) ^ (j + q - k)            (α_2) ^ (j + q - k)           ...          (α_(q - 2)) ^ (j + q - k)        0   1


- The matrix has q - k + 1 rows corresponding to the code's parity symbols.
- Any q - k + 1 columns form a Vandermonde matrix (non-singular).
- This ensures correction of up to q - k + 1 symbol erasures in a codeword.
- Permutation: We can re-arrange the columns of the `HField` matrix in any desired order.
- Parity symbols (s) deletion: Any set of `s` symbols within a codeword can be designated as parity symbols and permanently removed. This important property leads to construction of Shortened MDS codes.

Shortened MDS Codes:

- Corresponding columns of the field parity-check matrix `HField` can be deleted to form a shortened [[2 ^(m) + 1 - s, k, 2 ^ (m + 1) - s - k]] MDS code.
- This is an important property of MDS codes, particularly for their practical realisation in the form of augmented, extended RS codes because it enables efficient implementation in applications such as incremental redundancy systems, and network coding.
- 3 - level quantization of the received channel bits meaning 3 symbols deleted

Cyclic Code Construction:

- Using the first q - 1 columns of the field parity-check matrix (HField), using j = 0, and setting α_0, α_1, α_2, ..., α_q - 1 to  α ^ 0, α ^ 1, α ^ 2, ..., α ^ q - 1 in the parity-check matrix are set equal to the powers of a primitive element α of the Galois Field GF(q), a cyclic code can be constructed for efficient encoding and decoding. The resulting matrix is represented by `HSeed`.

- `HSeed` Matrix element expansion: 
    1. Row expansion: Each row of in the `HField` matrix is replaced with an `m`-by-`m` Field matrix defined over the base field GF (2 ^ m). This expansion is repesented by `HFieldExpanded`.
    2. Column expansion: The elements in each column of `HFieldExpanded` matrix are converted to binary representations by substituting powers of a primitive element (`α`) in the Galois Field GF(2 ^ m) with their corresponding m-tuples over the Boolean/Binary Field GF(2).
"""
function parity_checks(rs::ReedSolomon)
    r = ceil(Int, log2(rs.n + 1))
    GF2ͬ, a = finite_field(2, r, "a")
    q = 2^r + 1 - 3 # 3-level quantization
    HField = Matrix{FqFieldElem}(undef, q - rs.k + 1, q)
    for j in 1: q
        HField[1, j] = a ^ 0
    end
    HTemp2 = Matrix{FqFieldElem}(undef, r, q)
    for i in 1: q - rs.k + 1
        HField[i, 1] = a ^ 0
    end
    for i in 2:q - rs.k + 1
        for j in 2: q
            HField[i, j] = (a ^ (j - 1)) ^ (i - 2)
        end
    end
    HSeed = vcat(HField[1:1, :], HField[3:end, :])
    HFieldExpanded = Matrix{FqFieldElem}(undef, r * rs.k, q)
    g = 1
    while g <= r * rs.k
        for i in 1:q - rs.k
            for p in 1:r
                HTemp2[p:p, :] = reshape(HSeed[i, :].*a ^ (p - 1) , 1, :)
            end
        HFieldExpanded[g:g + r - 1, :] .=  HTemp2
        g = g + r
        end
    end
    H = Matrix{Bool}(undef, r * rs.k, r * q)
    for i in 1:r * rs.k
        for j in 1:q
            col_start = (j - 1) * r + 1
            col_end = col_start + r - 1
            t_tuple = Bool[]
            for k in 0:r - 1
                t_tuple = [t_tuple; !is_zero(coeff(HFieldExpanded[i, j], k))]
            end 
            H[i, col_start:col_end] .=  vec(t_tuple)
        end
    end
    return H 
end
