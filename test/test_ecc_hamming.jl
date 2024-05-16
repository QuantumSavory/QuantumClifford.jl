using Test
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, Hamming, Gilbert_Varshamov_bound, Hamming_bound, Griesmer_bound
using Nemo: matrix, GF

function designed_distance(matrix)
    for row in eachrow(matrix)
        count = sum(row)
        if count >= 3 # The minimum Hamming distance of the code is 3.
            return true
        end
    end
    return false
end

@testset "Testing Hamming codes properties" begin
    r_cases = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    for r in r_cases
        n = 2 ^ r - 1
        k = 2 ^ r - 1 - r
        H = parity_checks(Hamming(r))
        @test designed_distance(parity_checks(Hamming(r))) == true
        mat = matrix(GF(2), parity_checks(Hamming(r)))
        computed_rank = rank(mat)
        @test computed_rank == n - k
        @test Hamming_bound(Hamming(r)) == true
        @test Gilbert_Varshamov_bound(Hamming(r)) == true
    end
    # Example taken from [huffman2010fundamentals](@cite).
    @test parity_checks(Hamming(3)) == [0  0  0  1  1  1  1;
                                        0  1  1  0  0  1  1;
                                        1  0  1  0  1  0  1]
end
