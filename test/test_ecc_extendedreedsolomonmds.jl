using Test
using ILog2
using LinearAlgebra
using Combinatorics
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ExtendedReedSolomonMDS, generator_polynomial
using Nemo: finite_field, GF, FqFieldElem, FqPolyRingElem, coeff, is_zero, degree, matrix
include("utils_test_ecc.jl")

""" 
- Employing `3-level` quantization of the channel bits and erasing entire symbols if any of their constituent bits are erased can improve the performance of RS codes. Shortened Maximum Distance Separable (MDS) codes have the following parameters - code length `(n)`, codesize `(k)` and minimum Hamming distance `(d)` represented by [[n, k, d]] as follows: [[2 ^ (m) + 1 - s, k, 2 ^ (m + 1) - s - k]]. Thus, the designed minimum distance `d` is 2 ^ (m + 1) - s - k. Refer to chapter: 07, section: 03, pages: 172 to 175 [tomlinson2017error](@cite).
- The designed distance for binary expanded parity check matrix remains same as symbol based parity check matrix. According to [macwilliams1977theory](@cite), changing the basis `j` can increase the designed distance `(dmin)` of the resulting binary code.
"""

@testset "Testing designed distance of ExtendedReedSolomonMDS codes" begin
    m_cases = [3, 4, 5, 6, 7, 8]
    for m in m_cases  
        for t in rand(1:m - 1, 2)
            s_symbols = 3 # Refer to chapter: 07, section: 03, pages: 172 to 175 [tomlinson2017error](@cite).
            k = (2 ^ m -  1 - 2 * t) * m
            d = 2 ^ (m + 1) - s_symbols - k
            @test check_designed_distance(parity_checks(ExtendedReedSolomonMDS(m, t)), m, t, d, 0, 0) == true
        end
    end
end

@testset "Testing ExtendedReedSolomonMDS codes's properties" begin
    m_cases = [3, 4, 5, 6, 7, 8]
    for m in m_cases
        for t in rand(1:m - 1, 2)
            mat = matrix(GF(2), parity_checks(ExtendedReedSolomonMDS(m, t)))
            computed_rank = rank(mat) 
            s_symbols = 3 # Refer to chapter: 07, section: 03, pages: 172 to 175 [tomlinson2017error](@cite).
            k = (2 ^ m -  1 - 2 * t) * m
            n = (2 ^ m + 1 - s_symbols) * m
            @test computed_rank == n - k            
        end
    end

    # Example `(H₁₅₀₋₇₅`) taken from Eq. 7.9 of pg. 175 of [tomlinson2017error](@cite).
    @test size(parity_checks(ExtendedReedSolomonMDS(5, 8))) == (75, 150)
    H = Matrix{Bool}(undef, 75, 150)
    H = parity_checks(ExtendedReedSolomonMDS(5, 8))
    example₁₋₁ =	[1  0  0  0  0  1  0  0  0  0  1  0  0  0  0;	
			 0  1  0  0  0  0  1  0  0  0  0  1  0  0  0;
			 0  0  1  0  0  0  0  1  0  0  0  0  1  0  0;
			 0  0  0  1  0  0  0  0  1  0  0  0  0  1  0;
			 0  0  0  0  1  0  0  0  0  1  0  0  0  0  1;
			 1  0  0  0  0  0  1  0  0  0  0  0  1  0  0;
			 0  1  0  0  0  0  0  1  0  0  0  0  0  1  0;
			 0  0  1  0  0  0  0  0  1  0  0  0  0  0  1;
			 0  0  0  1  0  0  0  0  0  1  1  0  1  0  0;
			 0  0  0  0  1  1  0  1  0  0  0  1  0  1  0;
			 1  0  0  0  0  0  0  1  0  0  0  0  0  0  1;
			 0  1  0  0  0  0  0  0  1  0  1  0  1  0  0;
			 0  0  1  0  0  0  0  0  0  1  0  1  0  1  0;
			 0  0  0  1  0  1  0  1  0  0  0  0  1  0  1;
			 0  0  0  0  1  0  1  0  1  0  1  0  1  1  0;
			 1  0  0  0  0  0  0  0  1  0  0  1  0  1  0;
			 0  1  0  0  0  0  0  0  0  1  0  0  1  0  1;
			 0  0  1  0  0  1  0  1  0  0  1  0  1  1  0;
			 0  0  0  1  0  0  1  0  1  0  0  1  0  1  1;
			 0  0  0  0  1  0  0  1  0  1  1  0  0  0  1]
    example₁₋₂ =	[1  0  0  0  0  1  0  1  1  1  0  1  1  0  1;
			 0  1  0  0  0  1  1  1  1  1  1  0  0  1  0;
			 0  0  1  0  0  1  1  0  1  1  0  1  0  0  1;
			 0  0  0  1  0  1  1  0  0  1  1  0  0  0  0;
			 0  0  0  0  1  1  1  0  0  0  0  1  0  0  0]
    example₁₋₃ =	[1  0  0  0  0;
			 0  1  0  0  0; 
			 0  0  1  0  0;
			 0  0  0  1  0;
			 0  0  0  0  1;
			 1  0  0  1  0;
			 0  1  0  0  1;
			 1  0  0  0  0;
			 0  1  0  0  0;
			 0  0  1  0  0;
			 1  1  0  1  0;
			 0  1  1  0  1;
			 1  0  0  1  0;
			 0  1  0  0  1;
			 1  0  0  0  0;
			 1  0  0  1  1;
			 1  1  1  0  1;
			 1  1  0  1  0;
			 0  1  1  0  1;
			 1  0  0  1  0]
    example₁₋₄ =	[0  0  0  1  0; 
			 0  0  0  0  1;
			 1  0  1  0  0;
			 0  1  0  1  0;
			 0  0  1  0  1]

    @test H[1:20, 1:15] == example₁₋₁
    @test H[71:75, 1:15] == example₁₋₂ 
    @test H[1:20, 146:150] == example₁₋₃ 
    @test H[71:75, 146:150] == example₁₋₄
end
