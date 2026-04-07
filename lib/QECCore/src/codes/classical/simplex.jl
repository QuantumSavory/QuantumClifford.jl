"""
    $TYPEDEF

The `[2ʳ - 1, r, 2ʳ⁻¹]` simplex code family, dual to the binary Hamming codes.

`C(r)` is the dual of `Hamming(r)`. Its codewords are the rows of the Hamming
parity check matrix, and every nonzero codeword has weight `2ʳ⁻¹`.

Used as the seed code in the SHYPS construction [malcolm2025computing](@cite).
ECC Zoo: [Simplex code family](https://errorcorrectionzoo.org/c/simplex).

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
    # build the Hamming parity check matrix -- column j is just the number j written in r-bit binary
    H = zeros(Int, r, n)
    for j in 1:n
        for i in 1:r
            H[i, j] = (j >> (r - i)) & 1
        end
    end
    # the null space of H gives us the simplex parity check matrix
    # (because the simplex code is dual to Hamming, so its PCM = generator matrix of Hamming)
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
