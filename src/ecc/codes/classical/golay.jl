"""
The family of classical Golay binary codes were discovered by Edouard Golay in his 1949 paper [golay1949notes](@cite), where he described the binary `[[23, 12, 7]]` Golay code.

There are two binary Golay codes:

1. Binary `[[23, 12, 7]]` Golay Code: The perfect code with code length `(n)` 23 and dimension `(k)` 12. Minimum distance is 7, implying it can detect or correct up to 7 errors. By puncturing in any of the coordinates of parity check Matrix `H` = `[[24, 12, 8]]`, we obtain a `[[23, 12, 7]]` golay code.

2. Extended Binary [[24, 12, 8]] Golay Code: Obtained by adding a parity check bit to [[23, 12, 7]]. Thus, the binary codes can be obtained from each other by puncturing (shortening) or extending. 

Parity Check Matrix `(H)`: `H` is defined as follows: `H₂₄ = [ A' | I₁₂ ]` where I₁₂ is the 12 x 12 identity matrix and `A` is a bordered reverse circulant matrix.

Construction method for `A` [huffman2010fundamentals](@cite): The columns of A are labeled by ∞, 0, 1, 2, ..., 10. The first row contains 0 in column ∞ and 1 elsewhere. To obtain the second row, a 1 is placed in column ∞ nd a 1 is placed in columns 0, 1, 3, 4, 5, and 9; these numbers are precisely the squares of the integers modulo 11. That is, 0² = 0, 1² ≡ 10² ≡ 1 (mod 11), 2² ≡ 2² ≡ 4 (mod 11), etc. The third row of A is obtained by putting a 1 in column ∞ and then shifting the components in the second row one place to the left and wrapping the entry in column 0 around to column 10. The fourth row is obtained from the third in the same manner, as are the remaining rows.

Properties: The rows of `H₂₄` have weights 8 and 12. `H₂₄` is self-dual, meaning its codewords are also valid when their components are reversed. All codewords of `H₂₄` have weights divisible by four (`d` = 4 or 8). 

Punctured Code: All punctured codes are equivalent. Adding an overall parity check to `H₂₃` recovers `H₂₄`.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/golay).
"""
abstract type ClassicalCode end

struct Golay <: ClassicalCode
    n::Int 
   
    function Golay(n)
        if !(n ∈ [23, 24])
            throw(ArgumentError("Invalid parameters: `n` must be either 24 or 23 to obtain a valid code."))
        end
        new(n)
    end
end

function _circshift_row_golay(row::Vector{Int}, shift::Int, n::Int)
    l = length(row)
    return [row[mod((i - shift - 1), l) + 1] for i in 1:l]
end

function _create_A₂₄_golay(n::Int)
    A = zeros(Int, n ÷ 2, n ÷ 2)
    # Define the squared values modulo 11.
    squares_mod₁₁ = [0, 1, 4, 9, 5, 3, 3, 5, 9, 4, 1]
    A[1, 2:end] .= 1
    A[2, 1] = 1
    for i in squares_mod₁₁
        A[2, i + 2] = 1
    end
    # Fill in the rest of the rows using the reverse circulant property.
    for i in 3:n ÷ 2
        A[i, 2:end] = _circshift_row_golay(A[i - 1, 2:end], -1, n ÷ 2)
        A[i, 1] = 1  # Always put a 1 in the first column
    end
    return A
end

function parity_checks(g::Golay)
    if g.n == 24 # Extented binary `[[24, 12, 8]]` Golay code.
        A₂₄ = _create_A₂₄_golay(24)
        I₁₂ = Diagonal(ones(Int, g.n ÷ 2))
        H₂₄ = hcat(I₁₂, A₂₄')
        return H₂₄
    else # Perfect binary `[[23, 12, 7]]` Golay code.
        A₂₄ = _create_A₂₄_golay(24)
        A₂₃ = A₂₄[:, 1:end - 1]
        I₁₂ = Diagonal(ones(Int, g.n ÷ 2))
        H₂₃ = hcat(I₁₂, A₂₃')
        return H₂₃
    end
end
