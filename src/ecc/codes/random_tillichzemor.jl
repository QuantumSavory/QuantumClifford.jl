"""
The **random Tillich Zémor code** is a quantum LDPC code constructed using the
hypergraph product of two classical seed **(n, m, r)-Structured LDPC** codes.
"""
function random_TillichZemor_code end

function _create_matrix_M_random(rng::AbstractRNG, m::Int, n::Int, r::Int)
    M = zeros(Int, m, n - m)
    for col in 1:(n - m)
        # Randomly select r distinct rows to place ones
        rows = randperm(rng, m)[1:r]
        M[rows, col] .= 1
    end
    # Ensure no row is all zeros
    for row in 1:m
        if all(M[row, :] .== 0)
            # Randomly select a column and set this row to 1
            col = rand(rng, 1:(n - m))
            M[row, col] = 1
        end
    end
    return M
end

function _construct_parity_check_matrix(rng::AbstractRNG, n::Int, m::Int, r::Int)
    (m ≥ r && (n - m)*r ≥ m) || throw(ArgumentError(("Conditions for the existence of `M` in `H = [C | M]` are not satisfied.")))
    C = QECCore._create_circulant_matrix(m)
    M = _create_matrix_M_random(rng, m, n, r)
    # The parity-check matrix H = [C | M]
    H = hcat(C, M)
    return H
end

function random_TillichZemor_code(rng::AbstractRNG, n::Int, m::Int, r::Int)
    H = _construct_parity_check_matrix(rng, n, m, r)
    hx, hz = hgp(H, H)
    Stabilizer(CSS(hx, hz))
end

function random_TillichZemor_code(n::Int, m::Int, r::Int)
    H = _construct_parity_check_matrix(GLOBAL_RNG, n, m, r)
    hx, hz = hgp(H, H)
    Stabilizer(CSS(hx, hz))
end
