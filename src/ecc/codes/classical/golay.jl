"""
The family of classical binary Golay codes were discovered by Edouard Golay
in his 1949 paper [golay1949notes](@cite), where he described the binary
`[23, 12, 7]` Golay code.

There are two binary Golay codes:

- Binary `[23, 12, 7]` Golay code: The perfect code with code length `n = 23`
and dimension `k = 12`. By puncturing in any of the coordinates of parity check
matrix `H` = `[24, 12, 8]`, we obtain a `[23, 12, 7]` Golay code.

- Extended Binary `[24, 12, 8]` Golay code: Obtained by adding a parity check bit
to `[23, 12, 7]`. The bordered reverse circulant matrix `(A)` of `[24, 12, 8]`
Golay code is self-dual, i.e., `A₂₄` is same as A₂₄'.

The parity check matrix is defined as follows: `H₂₄ = [I₁₂ | A']` where `I₁₂` is the
`12 × 12` identity matrix and `A` is a bordered reverse circulant matrix. Puncturing
and then extending any column in​ with an overall parity check `H₂₃` reconstructs
the original parity check matrix `H₂₄`. Thus, all punctured codes are equivalent.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/golay).
"""
struct Golay <: ClassicalCode
    n::Int 
   
    function Golay(n)
        if !(n in (23, 24))
            throw(ArgumentError("Invalid parameters: `n` must be either 24 or 23 to obtain a valid code."))
        end
        new(n)
    end
end

function _circshift_row_golay(row::Vector{Int}, shift::Int, n::Int)
    l = length(row)
    return [row[mod((i - shift - 1), l) + 1] for i in 1:l]
end

# bordered reverse circulant matrix (see section 1.9.1, pg. 30-33) of [huffman2010fundamentals](@cite).
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
        A[i, 1] = 1 
    end
    return A
end

function generator(g::Golay)
    if g.n == 24 
        A₂₄ = _create_A₂₄_golay(24)
        I₁₂ = LinearAlgebra.Diagonal(ones(Int, g.n ÷ 2))
        G₂₄ = hcat(I₁₂, (A₂₄)')
        return G₂₄
    else 
        A₂₄ = _create_A₂₄_golay(24)
        A₂₃ = A₂₄[:, 1:end - 1]
        I₁₂ = LinearAlgebra.Diagonal(ones(Int, g.n ÷ 2))
        G₂₃ = hcat(I₁₂, (A₂₃)')
        return G₂₃
    end
end

function parity_checks(g::Golay)
    if g.n == 24
        A₂₄ = _create_A₂₄_golay(24)
        I₁₂ = LinearAlgebra.Diagonal(ones(Int, g.n ÷ 2))
        H₂₄ = hcat((A₂₄)', I₁₂)
        return H₂₄
    else 
        A₂₄ = _create_A₂₄_golay(24)
        A₂₃ = A₂₄[:, 1:end - 1]
        I₁₂ = LinearAlgebra.Diagonal(ones(Int, g.n ÷ 2))
        H₂₃ = hcat((A₂₃)', I₁₂)
        return H₂₃
    end
end

code_n(g::Golay) = g.n

code_k(g::Golay) = 12

distance(g::Golay) = code_n(g::Golay) - code_k(g::Golay) - 4
