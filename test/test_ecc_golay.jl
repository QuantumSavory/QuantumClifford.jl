using Test
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, Golay
using Nemo: matrix, GF

"""
- Theorem: Let `C` be a binary linear code. If `C` is self-orthogonal and has a generator matrix `G` where each row has weight divisible by four, then every codeword of `C` has weight divisible by four. 
- `H₂₄` is self-dual because its generator matrix has all rows with weight divisible by four. By above theorem, all codewords of `H₂₄` must have weights divisible by four. Refer to pages 30 to 33 of Ch1 of Fundamentals of Error Correcting Codes by Huffman, Cary and Pless, Vera.
""" 
function code_weight_property(matrix)
    for row in eachrow(matrix)
        count = sum(row)
        if count % 4 == 0
            return true
        end
    end
    return false
end

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
    @test code_weight_property(H) == true
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
