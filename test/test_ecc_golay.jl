using Test
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, Golay
using Nemo: matrix, GF

@testset "Testing binary Golay codes properties" begin
    test_cases = [(24, 12), (23, 12)]
    for (n, k) in test_cases
        H = parity_checks(Golay(n))
        mat = matrix(GF(2), parity_checks(Golay(n)))
        computed_rank = rank(mat)
        @test computed_rank == n - k
    end
    # [[24, 12, 8]] binary Golay code is a self-dual code [huffman2010fundamentals](@cite).
    H = parity_checks(Golay(24))
    @test H[:, (12 + 1):end]  == H[:, (12 + 1):end]'
    # Example taken from [huffman2010fundamentals](@cite).
    @test parity_checks(Golay(24)) == [1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1;
                                       0  1  0  0  0  0  0  0  0  0  0  0  1  1  1  0  1  1  1  0  0  0  1  0;
                                       0  0  1  0  0  0  0  0  0  0  0  0  1  1  0  1  1  1  0  0  0  1  0  1;
                                       0  0  0  1  0  0  0  0  0  0  0  0  1  0  1  1  1  0  0  0  1  0  1  1;
                                       0  0  0  0  1  0  0  0  0  0  0  0  1  1  1  1  0  0  0  1  0  1  1  0;
                                       0  0  0  0  0  1  0  0  0  0  0  0  1  1  1  0  0  0  1  0  1  1  0  1;
                                       0  0  0  0  0  0  1  0  0  0  0  0  1  1  0  0  0  1  0  1  1  0  1  1;
                                       0  0  0  0  0  0  0  1  0  0  0  0  1  0  0  0  1  0  1  1  0  1  1  1;
                                       0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  1  0  1  1  0  1  1  1  0;
                                       0  0  0  0  0  0  0  0  0  1  0  0  1  0  1  0  1  1  0  1  1  1  0  0;
                                       0  0  0  0  0  0  0  0  0  0  1  0  1  1  0  1  1  0  1  1  1  0  0  0;
                                       0  0  0  0  0  0  0  0  0  0  0  1  1  0  1  1  0  1  1  1  0  0  0  1]

    # Example taken from https://mathworld.wolfram.com/GolayCode.html.
    @test parity_checks(Golay(23)) == [1  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1;
                                       0  1  0  0  0  0  0  0  0  0  0  1  1  1  0  1  1  1  0  0  0  1  0;
                                       0  0  1  0  0  0  0  0  0  0  0  1  1  0  1  1  1  0  0  0  1  0  1;
                                       0  0  0  1  0  0  0  0  0  0  0  1  0  1  1  1  0  0  0  1  0  1  1;
                                       0  0  0  0  1  0  0  0  0  0  0  1  1  1  1  0  0  0  1  0  1  1  0;
                                       0  0  0  0  0  1  0  0  0  0  0  1  1  1  0  0  0  1  0  1  1  0  1;
                                       0  0  0  0  0  0  1  0  0  0  0  1  1  0  0  0  1  0  1  1  0  1  1;
                                       0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  1  0  1  1  0  1  1  1;
                                       0  0  0  0  0  0  0  0  1  0  0  1  0  0  1  0  1  1  0  1  1  1  0;
                                       0  0  0  0  0  0  0  0  0  1  0  1  0  1  0  1  1  0  1  1  1  0  0;
                                       0  0  0  0  0  0  0  0  0  0  1  1  1  0  1  1  0  1  1  1  0  0  0]
end
