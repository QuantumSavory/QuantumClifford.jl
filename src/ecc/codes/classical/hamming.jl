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
struct Hamming <: AbstractCECC
    r::Int
    function Hamming(r)
        if r < 2
            throw(ArgumentError("Invalid parameters: `r` must be ≥ 2 to obtain a valid code."))
        end
        new(r)
    end
end

function parity_matrix(h::Hamming)
    n = 2^h.r - 1
    max_elem = n * h.r
    rows = Vector{Int}(undef, max_elem)
    cols = Vector{Int}(undef, max_elem)
    vals = Vector{Int}(undef, max_elem)
    idx = 1
    @inbounds for j in 1:n
        mask = 1 << (h.r - 1)
        @simd for i in 1:h.r
            if j & mask != 0
                rows[idx] = i
                cols[idx] = j
                vals[idx] = 1
                idx += 1
            end
            mask >>= 1
        end
    end
    rows = rows[1:idx-1]
    cols = cols[1:idx-1]
    vals = vals[1:idx-1]
    return sparse(rows, cols, vals, h.r, n)
end

code_n(h::Hamming) = 2 ^ h.r - 1

code_k(h::Hamming) = 2 ^ h.r - 1 - h.r

distance(h::Hamming) = 3
