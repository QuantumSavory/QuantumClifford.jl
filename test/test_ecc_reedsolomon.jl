using Test
using Nemo
using Combinatorics
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedSolomon

test_cases = [(5, 3), (8, 2), (6, 4), (10, 5), (7, 2)]
@testset "# Test for the defining property of a Reed-Solomon codes" begin
     for (i, (m, e)) in enumerate(test_cases)
         G = generator_matrix(ReedSolomon(m, e))
         @test rank(G) == m 
     end
end
