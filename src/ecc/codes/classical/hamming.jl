using SparseArrays

abstract type ClassicalCode end

"""
The family of `[2ʳ - 1, 2ʳ - 1 - r, 3]` Hamming binary codes were discovered
by Richard W. Hamming in his 1950 paper [hamming1950error](@cite) as a way of
automatically correcting errors introduced by punched card readers. In his
original paper, Hamming elaborated his general idea, but specifically focused
on the `Hamming(7, 4)` code which adds three parity bits to four bits of data.

The Hamming matrix `H` is an `r × (2ʳ - 1)` binary matrix, where each column
corresponds to the binary representation of the integers from 1 to 2ʳ - 1, with
`r ≥ 2`. This matrix serves as the parity-check matrix for a binary Hamming code
with parameters `[2ʳ − 1, 2ʳ − 1 − r, 3]` The minimum Hamming distance of this
code is 3, as detailed in [huffman2010fundamentals](@cite).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/hamming).
"""
struct Hamming <: ClassicalCode
    r::Int
    function Hamming(r)
        if r < 2
            throw(ArgumentError("Invalid parameters: `r` must be ≥ 2 to obtain a valid code."))
        end
        new(r)
    end
end

function parity_checks(h::Hamming)
    n = 2^h.r - 1 # Number of columns in H
    max_elements = n * h.r # Max non-zero entries in H
    # Pre-allocate arrays for sparse matrix indices and values
    rows = Vector{Int}(undef, max_elements)
    cols = Vector{Int}(undef, max_elements)
    vals = Vector{Int}(undef, max_elements)
    idx = 1 # Tracks position in arrays
    @inbounds for j in 1:n
        mask = 1 << (h.r - 1) # Initialize mask for MSB
        @simd for i in 1:h.r
            if j & mask != 0 # Check if the i-th bit is 1
                rows[idx] = i
                cols[idx] = j
                vals[idx] = 1
                idx += 1
            end
            mask >>= 1 # Shift mask to next bit
        end
    end
    # Resize arrays to actual number of non-zero elements
    rows = rows[1:idx-1]
    cols = cols[1:idx-1]
    vals = vals[1:idx-1]
    return sparse(rows, cols, vals, h.r, n)
end

code_n(h::Hamming) = 2 ^ h.r - 1

code_k(h::Hamming) = 2 ^ h.r - 1 - h.r

distance(h::Hamming) = 3
