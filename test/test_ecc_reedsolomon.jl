using Test
using ILog2
using LinearAlgebra
using Combinatorics
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedSolomon, generator_polynomial
using Nemo: finite_field, GF, FqFieldElem, FqPolyRingElem, coeff, is_zero, degree, matrix
include("utils_test_ecc.jl")

""" 
- Employing `3-level` quantization of the channel bits and erasing entire symbols if any of their constituent bits are erased can improve the performance of RS codes. Shortened Maximum Distance Separable (MDS) codes have the following parameters - code length `(n)`, codesize `(k)` and minimum Hamming distance `(d)` represented by [[n, k, d]] as follows: [[2 ^ (m) + 1 - s, k, 2 ^ (m + 1) - s - k]]. Thus, the designed minimum distance `d` is 2 ^ (m + 1) - s - k. Refer to chapter: 07, section: 03, pages: 172 to 175 [tomlinson2017error](@cite).
- The designed distance for binary expanded parity check matrix remains same as symbol based parity check matrix. According to [macwilliams1977theory](@cite), changing the basis `j` can increase the designed distance `(dmin)` of the resulting binary code.
"""

@testset "Testing designed distance of Shortened and Maximum Distance Separable (MDS) Reed-Solomon codes" begin
    m_cases = [3, 4, 5, 6, 7, 8]
    for m in m_cases  
        for t in rand(1:m - 1, 2)
            # The variable s_symbols is set to 3 because the last three columns of the field parity check matrix `HF`​ are erased. By using the first `x − 1` columns of `HF`​, assigning `j = 0` and setting `α₀, α₁, α₂, ..., αₓ ₋ ₁` to `α⁰, α¹, α², ..., αˣ ⁻ ¹`, where `α` is a primitive element of `GF(x)`. This approach allows for the construction of a cyclic code, which offers advantages in the implementation of encoding and decoding.  
            s_symbols = 3 
            k = (2 ^ m -  1 - 2 * t) * m
            d = 2 ^ (m + 1) - s_symbols - k
            @test check_designed_distance(parity_checks(ReedSolomon(m, t)), m, t, d, 0, 0) == true
        end
    end
end

@testset "Testing Shortened and Maximum Distance Separable (MDS) Reed Solomon codes's properties" begin
    m_cases = [3, 4, 5, 6, 7, 8]
    for m in m_cases
        for t in rand(1:m - 1, 2)
            mat = matrix(GF(2), parity_checks(ReedSolomon(m, t)))
            computed_rank = rank(mat) 
            s_symbols = 3 
            k = (2 ^ m -  1 - 2 * t) * m
            n = (2 ^ m + 1 - s_symbols) * m
            @test computed_rank == n - k            
            n_gx = 2 ^ m - 1
            GF2ʳ, a = finite_field(2, m, "a")
            GF2x, x = GF2ʳ[:x]
            # Reed-Solomon code is cyclic as its generator polynomial, `g(x)` divides `xⁿ - 1`, so `mod (xⁿ - 1, g(x))` = 0.
            @test mod(x ^ n_gx - 1, generator_polynomial(ReedSolomon(m, t))) == 0
        end
    end

    # Examples taken from pg. 18 of http://hscc.cs.nthu.edu.tw/~sheujp/lecture_note/rs.pdf.
    GF2ʳ, a = finite_field(2, 3, "a")
    P, x = GF2ʳ[:x]
    @test generator_polynomial(ReedSolomon(3, 2)) == x ^ 4 + (a + 1) * x ^ 3 + x ^ 2 + a * x + a + 1
    
    # Example taken from https://www.youtube.com/watch?v=dpxD8gwgbOc.
    GF2ʳ, a = finite_field(2, 4, "a")
    P, x = GF2ʳ[:x]
    @test generator_polynomial(ReedSolomon(4, 2)) == x ^ 4 + a ^ 13 * x ^ 3 + a ^ 6 * x ^ 2 + a ^ 3 * x + a ^ 10

    # Example `(H₁₅₀₋₇₅`) taken from Eq. 7.9 of pg. 175 of [tomlinson2017error](@cite).
    @test size(parity_checks(ReedSolomon(5, 8))) == (75, 150)
    H = Matrix{Bool}(undef, 75, 150)
    H = parity_checks(ReedSolomon(5, 8))
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
