using Test
using Combinatorics
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, QuantumReedMuller

function split(H::Matrix{Bool})
    n = size(H, 2)
    m = div(n, 2)
    Hx = H[:, 1:m]
    Hz = H[:, m+1:end]
    return Hx, Hz
end

function min_distance(Hx::Matrix{Bool}, Hz::Matrix{Bool})
    w = typemax(Int)
    for (x, z) in zip(eachrow(Hx), eachrow(Hz))
        r = map(|, x, z)
        w = min(w, sum(r))
    end
    return w 
end

@testset "Test QRM(t, r) Matrix Minimum Distance" begin
    # Prepare tests for different parameters
    for (t, r) in [(1, 3), (2, 4)]
        # Generate the parity check matrix
        stab = parity_checks(QuantumReedMuller(t, r))
        H = stab_to_gf2(stab)
        Hx, Hz = split(H)
        computed_distance = min_distance(Hx, Hz)
        expected_distance = 2^t + 2^(t - 1)
        @test computed_distance == expected_distance
    end
end