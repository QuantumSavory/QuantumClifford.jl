using Test
using ILog2
using LinearAlgebra
using Combinatorics
using Nemo: finite_field, GF, FqFieldElem, FqPolyRingElem, coeff, is_zero, degree, matrix
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedSolomon, generator_polynomial

""" 
- Employing `3-level` quantization of the channel bits and erasing entire symbols if any of their constituent bits are erased can improve the performance of RS codes. Shortened Maximum Distance Separable (MDS) codes have the following parameters - code length `(n)`, codesize `(k)` and minimum Hamming distance `(d)` represented by [[n, k, d]] as follows: [[2 ^ (m) + 1 - s, k, 2 ^ (m + 1) - s - k]]. Thus, the designed minimum distance `d` is 2 ^ (m + 1) - s - k. Refer to chapter: 07, section: 03, pages: 172 to 175 [tomlinson2017error](@cite).
- The designed distance for binary expanded parity check matrix remains same as symbol based parity check matrix. According to [macwilliams1977theory](@cite), changing the basis `j` can increase the designed distance `(dmin)` of the resulting binary code.
"""
function check_designed_distance(matrix, m, t)
    k = 2 ^ m -  1 - 2 * t
    s_symbols = 3
    n_cols = size(matrix, 2)
    for num_cols in 1:2 ^ (m + 1) - s_symbols - k
        for i in 1:n_cols - num_cols + 1
            combo = matrix[:, i:(i + num_cols - 1)]
            sum_cols = sum(combo, dims = 2)
            if all(sum_cols .== 0)
                return false  # Minimum distance is not greater than `2 ^ (m + 1) - s_symbols - k`.
            end
        end
    end
    return true  # Minimum distance is at least `2 ^ (m + 1) - s_symbols - k`.
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
            d = 2 ^ (m + 1) - s_symbols - k
            @test computed_rank == n - k
            @test check_designed_distance(parity_checks(ReedSolomon(m, t)), m, t) == true
            # Reed-Solomon codes exactly meet the [Singleton Bound](https://en.wikipedia.org/wiki/Singleton_bound).
            @test d <= n - k + 1 
            # Reed-Solomon code is cyclic as its generator polynomial, `g(x)` divides `xⁿ - 1`, so `mod (xⁿ - 1, g(x))` = 0.
            n_gx = 2 ^ m - 1
            GF2ʳ, a = finite_field(2, m, "a")
            GF2x, x = GF2ʳ[:x]
            mod(x ^ n_gx - 1, generator_polynomial(ReedSolomon(m, t))) == 0
        end
    end

    # RS(7, 3), RS(15, 9), RS(255, 223), RS(160, 128), RS(255, 251), (255, 239) and (255, 249) codes. Examples taken from https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction, https://www.cs.cmu.edu/~guyb/realworld/reedsolomon/reed_solomon_codes.html, http://www.chencode.cn/lecture/Information_Theory_and_Coding/Information%20Theory%20and%20Coding-CH7.pdf, https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=91e1d6d27311780b0a8c34a41793fa85f3947af1.
    test_cases = [(7, 3), (15, 9), (225, 223), (160, 128), (255, 251), (255, 239), (255, 249)]
    for (n, k) in test_cases
        m = ilog2(n + 1)
        t = div(n - k, 2)
        # Using fixed generator polynomial construction scheme for defining generator polynomial, `g(x)`, of RS codes, `degree(g(x))` == 2 * t == n - k. 
        @test degree(generator_polynomial(ReedSolomon(m, t))) == 2 * t == n - k
    end

    # Examples taken from [http://hscc.cs.nthu.edu.tw/~sheujp/lecture_note/rs.pdf].
    GF2ʳ, a = finite_field(2, 3, "a")
    P, x = GF2ʳ[:x]
    @test generator_polynomial(ReedSolomon(3, 2)) == x ^ 4 + (a + 1) * x ^ 3 + x ^ 2 + a * x + a + 1
    
    # Example taken from https://www.youtube.com/watch?v=dpxD8gwgbOc
    GF2ʳ, a = finite_field(2, 4, "a")
    P, x = GF2ʳ[:x]
    @test generator_polynomial(ReedSolomon(4, 2)) == x ^ 4 + a ^ 13 * x ^ 3 + a ^ 6 * x ^ 2 + a ^ 3 * x + a ^ 10

    # Example `(H₁₅₀₋₇₅`) taken from Eq. 7.9 of pg. 175 of [tomlinson2017error](@cite).
    @test size(parity_checks(ReedSolomon(5, 8))) == (75, 150)
    H = Matrix{Bool}(undef, 75, 150)
    H = parity_checks(ReedSolomon(5, 8))
    example₁₋₁ = [1  0  0  0  0  1  0  0  0  0  1  0  0  0  0;	
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
    example₁₋₂ = [1  0  0  0  0  1  0  1  1  1  0  1  1  0  1;
			 0  1  0  0  0  1  1  1  1  1  1  0  0  1  0;
			 0  0  1  0  0  1  1  0  1  1  0  1  0  0  1;
			 0  0  0  1  0  1  1  0  0  1  1  0  0  0  0;
			 0  0  0  0  1  1  1  0  0  0  0  1  0  0  0]
    example₁₋₃ = [1  0  0  0  0;
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
    example₁₋₄ = [0  0  0  1  0; 
			 0  0  0  0  1;
			 1  0  1  0  0;
			 0  1  0  1  0;
			 0  0  1  0  1]

    @test H[1:20, 1:15] == example₁₋₁
    @test H[71:75, 1:15] == example₁₋₂ 
    @test H[1:20, 146:150] == example₁₋₃ 
    @test H[71:75, 146:150] == example₁₋₄
end
