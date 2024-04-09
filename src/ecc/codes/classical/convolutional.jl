"""The family of Convolutional codes, as discovered by Elias in his 1955 paper [caldwell1955processing](@cite).

You might be interested in consulting [lieb2021convolutional](@cite), [viterbi1967error](@cite), and [forney1970convolutional](@cite) and [forney1974convolutional](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/convolutional)
"""
struct Convolutional <: ClassicalCode
    n::Int
    k::Int
    q::Int

    function Convolutional(n, k, q)
        if n < 0 || k < 0 && k <= n|| q < 0 || n > 20 || k > 20 || q > 20
            throw(ArgumentError("Invalid parameters: n, k, q must be non-negative  and <= 20 in order to obtain a valid code and to remain tractable"))
        end
        new(n, k, q)
    end
end

# Formula: C = im(G) = {uG∣ u ∈ Fₓ [D]ᵏ} 
function _create_generator_matrix(P, x, k, n)
    G = Matrix{FqPolyRingElem}(undef, k, n)
    for i in 1:k
        for j in 1:n
            if j >= i
                G[i, j] = x^(j - i)
            else
                G[i, j] = zero(P)
            end
        end
    end
    return G
end

# Formula: C = ker⁡(H) = {v ∣ v ∈ Fₓ[D]ⁿ, vHᵗ = 0}
function _create_parity_check_matrix(P, k, n)
    H = Matrix{FqPolyRingElem}(undef, n - k, n)
    for i in 1:(n - k)
        for j in 1:n
            if j <= i
                H[i, j] = one(P)
            else
                H[i, j] = zero(P)
            end
        end
    end
    return H
end

function generator_polynomial(rs::Convolutional)
    r = ceil(Int, log2(rs.n + 1))
    GF2q, a = finite_field(rs.q, r , "a")
    P, x = GF2q[:x]
    G = _create_generator_matrix(P, x, rs.k, rs.n)
    H = _create_parity_check_matrix(P, rs.k, rs.n)
    return G, H
end
