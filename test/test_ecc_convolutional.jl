using Test
using Nemo
using Combinatorics
using AbstractAlgebra
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC

# Formula: C = im(G) = {uG∣ u ∈ Fₓ [D]ᵏ} 
function create_generator_matrix(P, x, k, n)
    G = Array{FqFieldElem}(undef, k, n)
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


# Formula: C = im(G) = {uG∣ u ∈ Fₓ [D]ᵏ} 
function create_generator_matrix(P, x, k, n)
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
function create_parity_check_matrix(P, k, n)
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


@testset "Test Convolutional Codes" begin
    k_values = [2, 3, 4] 
    n_values = [6, 8, 10] 
    q = 2
    for k in k_values, n in n_values
        r = ceil(Int, log2(n + 1))
        GF2, a = finite_field(q, r , "a")
        P, x = GF2[:x]
        G = create_generator_matrix(P, x, k, n)
        H = create_parity_check_matrix(P, k, n)
        @test size(G) == (k, n)  
        @test size(H) == (n - k, n)
    end
end
