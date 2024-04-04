using Test
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedSolomon

@testset "Test ReedSolomon(m, e) Matrix Rank" begin
    for m in 5:20
        for e in 2:10
            G = generator_matrix(ReedSolomon(m, e))
            @test rank(G) == m
        end
    end
end
