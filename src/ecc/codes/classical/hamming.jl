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
        if r <= 2
            throw(ArgumentError("Invalid parameters: `r` must be ≥ 2 to obtain a valid code."))
        end
        new(r)
    end
end

using SparseArrays

function parity_checks(h::Hamming)
    n = 2 ^ h.r - 1
    rows = Int[]
    cols = Int[]
    vals = Int[]
    for j in 1:n
        columnsⱼ = bitstring(j)[end - h.r + 1:end]
        for i in 1:h.r
            if parse(Int, columnsⱼ[i]) == 1
                push!(rows, i)
                push!(cols, j)
                push!(vals, 1)
            end
        end
    end
    H = sparse(rows, cols, vals, h.r, n)
    return H
end

code_n(h::Hamming) = 2 ^ h.r - 1

code_k(h::Hamming) = 2 ^ h.r - 1 - h.r

distance(h::Hamming) = 3
