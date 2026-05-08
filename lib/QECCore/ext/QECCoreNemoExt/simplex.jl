"""
    $TYPEDEF

The `[2ʳ - 1, r, 2ʳ⁻¹]` simplex code family, dual to the binary Hamming codes.

`C(r)` is the dual of `Hamming(r)`. Its codewords are the rows of the Hamming
parity check matrix, and every nonzero codeword has weight `2ʳ⁻¹`.

Used as the seed code in the SHYPS construction [malcolm2025computing](@cite).
ECC Zoo: [Simplex code family](https://errorcorrectionzoo.org/c/simplex).

### Fields
    \$TYPEDFIELDS
"""
struct Simplex <: AbstractCECC
    """Parameter `r` (must be ≥ 2). Determines code length `n = 2ʳ - 1`."""
    r::Int
    function Simplex(r)
        if r < 2
            throw(ArgumentError("Invalid parameters: `r` must be ≥ 2 to obtain a valid code."))
        end
        new(r)
    end
end

"""
    _dual(H)

Compute the dual code of a binary parity check matrix `H` using Nemo's nullspace.
Returns the parity check matrix of the dual code as a transposed nullspace matrix.

This is a general-purpose utility: for any binary matrix `H`, the dual code's
generator matrix `G` satisfies `H * Gᵀ = 0` over GF(2).
"""
function _dual(H)
    r, c = size(H)
    R, _ = Nemo.residue_ring(Nemo.ZZ, 2)
    H_nemo = Nemo.matrix_space(R, r, c)(R.(H))
    null = Nemo.nullspace(H_nemo)[2]
    @assert all(iszero, H_nemo * null)
    return Nemo.transpose(null)
end

function QECCore.parity_matrix(c::Simplex)
    r = c.r
    n = 2^r - 1
    # Building the Hamming parity check matrix -- column j is the number j in r-bit binary
    H_hamming = zeros(Int, r, n)
    for j in 1:n
        for i in 1:r
            H_hamming[i, j] = (j >> (r - i)) & 1
        end
    end
    # The dual of the Hamming code is the Simplex code
    dual_mat = _dual(H_hamming)
    # Converting Nemo matrix back to Int matrix
    nr, nc = size(dual_mat)
    result = zeros(Int, nr, nc)
    F1 = Nemo.residue_ring(Nemo.ZZ, 2)[1](1)
    for i in 1:nr, j in 1:nc
        result[i, j] = dual_mat[i, j] == F1 ? 1 : 0
    end
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

QECCore.code_n(c::Simplex) = 2^c.r - 1
QECCore.distance(c::Simplex) = 2^(c.r - 1)
