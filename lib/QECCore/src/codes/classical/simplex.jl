"""
    $TYPEDEF

The family of `[2ʳ - 1, r, 2ʳ⁻¹]` simplex codes, dual to the binary Hamming codes.

The simplex code `C(r)` is the dual of the Hamming code `Hamming(r)`. Its codewords
are the rows of the Hamming parity check matrix, and every nonzero codeword has
weight exactly `2ʳ⁻¹`. The code has parameters `[2ʳ − 1, r, 2ʳ⁻¹]`.

The simplex code is used as the seed code for the Subsystem Hypergraph Product
Simplex (SHYPS) quantum code construction [malcolm2025computing](@cite).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/simplex).

### Fields
    $TYPEDFIELDS
"""
struct Simplex <: AbstractCECC
    r::Int
    function Simplex(r)
        if r < 2
            throw(ArgumentError("Invalid parameters: `r` must be ≥ 2 to obtain a valid code."))
        end
        new(r)
    end
end

"""
    _gf2_nullspace(M::AbstractMatrix{<:Integer})

Compute the null space of a binary matrix `M` over GF(2) using Gaussian elimination.
Returns a matrix whose rows form a basis for the null space.
"""
function _gf2_nullspace(M::AbstractMatrix{<:Integer})
    m, n = size(M)
    A = mod.(Matrix(M), 2)
    pivot_cols = Int[]
    pivot_row = 1
    for j in 1:n
        pivot = 0
        for i in pivot_row:m
            if A[i, j] == 1
                pivot = i
                break
            end
        end
        pivot == 0 && continue
        if pivot != pivot_row
            A[pivot_row, :], A[pivot, :] = A[pivot, :], A[pivot_row, :]
        end
        push!(pivot_cols, j)
        for i in 1:m
            if i != pivot_row && A[i, j] == 1
                @. A[i, :] = mod(A[i, :] + A[pivot_row, :], 2)
            end
        end
        pivot_row += 1
    end
    free_cols = setdiff(1:n, pivot_cols)
    nullity = length(free_cols)
    if nullity == 0
        return zeros(Int, 0, n)
    end
    N = zeros(Int, nullity, n)
    for (k, fc) in enumerate(free_cols)
        N[k, fc] = 1
        for (pr, pc) in enumerate(pivot_cols)
            N[k, pc] = A[pr, fc]
        end
    end
    return N
end

function parity_matrix(c::Simplex)
    r = c.r
    n = 2^r - 1
    # Build Hamming PCM H (r × n): columns are binary representations of 1..n
    H = zeros(Int, r, n)
    for j in 1:n
        for i in 1:r
            H[i, j] = (j >> (r - i)) & 1
        end
    end
    # Null space of Hamming PCM = generator matrix of Hamming code = PCM of Simplex code
    result = _gf2_nullspace(H)
    rows_idx = Int[]
    cols_idx = Int[]
    for i in 1:size(result, 1), j in 1:size(result, 2)
        if result[i, j] == 1
            push!(rows_idx, i)
            push!(cols_idx, j)
        end
    end
    return sparse(rows_idx, cols_idx, ones(Int, length(rows_idx)), size(result, 1), n)
end

code_n(c::Simplex) = 2^c.r - 1

code_k(c::Simplex) = c.r

distance(c::Simplex) = 2^(c.r - 1)
