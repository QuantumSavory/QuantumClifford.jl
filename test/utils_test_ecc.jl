"""The designed distance of a classical linear code with parameters `[[n, k, d]]`, where `n` represents the code length, `k` denotes the code dimension, and `d` signifies the minimum Hamming distance. For polynomial codes, `t` indicates the degree of the generator polynomial, and `m` represents the extension degree for the finite Galois field `GF(2ᵐ)`."""
function check_designed_distance(matrix, m, t, d, n, k)
    n_cols = size(matrix, 2)
    for num_cols in 1:d
        for i in 1:n_cols - num_cols + 1
            combo = matrix[:, i:(i + num_cols - 1)]
            sum_cols = sum(combo, dims = 2)
            if all(sum_cols .== 0)
                return false 
            end
        end
    end
    return true 
end
