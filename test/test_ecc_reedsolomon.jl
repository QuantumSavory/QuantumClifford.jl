using Test
using Combinatorics
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedSolomon, generator_polynomial
using Nemo: finite_field, GF, FqFieldElem, FqPolyRingElem, coeff, is_zero, degree, matrix

@testset "Testing ReedSolomon codes's properties" begin
    m_cases = [3, 4, 5, 6, 7, 8]
    for m in m_cases
        for t in rand(1:m - 1, 2)   
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
end
