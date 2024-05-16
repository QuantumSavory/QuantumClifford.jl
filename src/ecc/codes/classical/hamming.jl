"""
The family of `[[2ʳ - 1, 2ʳ - 1 - r, 3]]` Hamming binary codes were discovered by Richard W. Hamming in his 1950 paper [hamming1950error](@cite) as a way of automatically correcting errors introduced by punched card readers. In his original paper, Hamming elaborated his general idea, but specifically focused on the Hamming(7, 4) code which adds three parity bits to four bits of data.

Let `n = 2ʳ - 1`, with `r ≥ 2`. Then the `r × (2ʳ - 1)` matrix `H` whose columns, in order, are the numbers `1, 2, . . . , 2ʳ - 1` written as binary numerals, is the parity check matrix of an `[n = 2ʳ − 1, k = n − r]` binary code [huffman2010fundamentals](@cite). Any rearrangement of columns of Hᵣ gives an equivalent code. The code is denoted by either `Hᵣ` or `H(2, r)`. The minimum Hamming distance of these code is 3.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/hamming).
"""

abstract type ClassicalCode end

struct Hamming <: ClassicalCode
    r::Int 
   
    function Hamming(r)
        if r <= 2
            throw(ArgumentError("Invalid parameters: `r` must be ≥ 2 to obtain a valid code."))
        end
        new(r)
    end
end

function parity_checks(h::Hamming)
    n = 2 ^ h.r - 1
    H = zeros(Int, h.r, n)
    for j in 1:n
        columnsⱼ = bitstring(j)[end - h.r + 1:end]
        for i in 1:h.r
            H[i, j] = parse(Int, columnsⱼ[i])
        end
    end
    return H
end

code_n(h::Hamming) = 2 ^ h.r - 1
code_k(h::Hamming) = 2 ^ h.r - 1 - h.r
distance(h::Hamming) = 3
