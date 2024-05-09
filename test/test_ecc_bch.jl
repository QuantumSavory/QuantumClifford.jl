using Test
using Nemo: ZZ, residue_ring, matrix, finite_field, GF, minpoly, coeff, lcm, FqPolyRingElem, FqFieldElem, is_zero, degree, defining_polynomial, is_irreducible 
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, BCH, generator_polynomial

"""
- To prove that `t`-bit error correcting BCH code indeed has minimum distance of at least `2 * t + 1`, it is shown that no `2 * t` or fewer columns of binary parity check matrix `H` sum to zero. A formal mathematical proof can be found on pages 202 and 203 of https://personal.oss.unist.hr/~mnizetic/ZASTITNO%20LINIJSKI%20KODIRANJE/SEMINARSKI%20RADOVI/CH06.pdf. 
- The paramater `2 * t + 1` is usually called the designed distance of the `t`-bit error correcting BCH code.
"""
function check_designed_distance(matrix, t)
    n_cols = size(matrix, 2)
    for num_cols in 1:2 * t
        for i in 1:(n_cols - num_cols + 1)
            combo = matrix[:, i:(i + num_cols - 1)]
            sum_cols = sum(combo, dims = 2)
            if all(sum_cols .== 0)
                return false  # Minimum distance is not greater than 2 * t
            end
        end
    end
    return true  # Minimum distance is at least 2 * t + 1
end

@testset "Testing properties of BCH codes" begin
    m_cases = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    for m in m_cases
        n = 2 ^ m - 1
        for t in 1:m
            H = parity_checks(BCH(m, t))
            @test check_designed_distance(H, t) == true
            # n - k == degree of generator polynomial, `g(x)` == rank of binary parity check matrix, `H`.
            mat = matrix(GF(2), parity_checks(BCH(m, t)))
            computed_rank = rank(mat)
            @test computed_rank == degree(generator_polynomial(BCH(m, t)))
            @test code_k(BCH(m, t)) == n - degree(generator_polynomial(BCH(m, t)))
            # BCH code is cyclic as its generator polynomial, `g(x)` divides `xⁿ - 1`, so `mod (xⁿ - 1, g(x))` = 0.
            gx = generator_polynomial(BCH(m, t))
            GF2x, x = GF(2)["x"] 
            @test mod(x ^ n - 1, gx) == 0
        end
    end

    #example taken from https://personal.oss.unist.hr/~mnizetic/ZASTITNO%20LINIJSKI%20KODIRANJE/SEMINARSKI%20RADOVI/CH06.pdf.
    @test parity_checks(BCH(4, 2))  ==    [1  0  0  0  1  0  0  1  1  0  1  0  1  1  1;
                                           0  1  0  0  1  1  0  1  0  1  1  1  1  0  0;
                                           0  0  1  0  0  1  1  0  1  0  1  1  1  1  0;
                                           0  0  0  1  0  0  1  1  0  1  0  1  1  1  1;
                                           1  0  0  0  1  1  0  0  0  1  1  0  0  0  1;
                                           0  0  0  1  1  0  0  0  1  1  0  0  0  1  1;
                                           0  0  1  0  1  0  0  1  0  1  0  0  1  0  1;
                                           0  1  1  1  1  0  1  1  1  1  0  1  1  1  1]
    GF2x, x = GF(2)["x"]
    GF2⁴, a = finite_field(2, 4, "a")
    GF2⁶, b = finite_field(2, 6, "b")
    #examples taken from https://web.ntpu.edu.tw/~yshan/BCH_code.pdf.
    @test defining_polynomial(GF2x, GF2⁴) ==  x ^ 4 + x + 1
    @test is_irreducible(defining_polynomial(GF2x, GF2⁴)) == true
    @test generator_polynomial(BCH(4, 2)) == x ^ 8 + x ^ 7 +  x ^ 6 +  x ^ 4 + 1
    @test generator_polynomial(BCH(4, 3)) == x ^ 10 + x ^ 8 +  x ^ 5 +  x ^ 4 +  x ^ 2 + x + 1
    @test defining_polynomial(GF2x, GF2⁶) == x ^ 6 + x ^ 4 + x ^ 3 + x + 1
    @test is_irreducible(defining_polynomial(GF2x, GF2⁶)) == true
    # GF(2ⁿ) can have multiple distinct primitive polynomials of the same degree due to existence of several irreducible polynomials of that degree, each generating the field through different roots.
    # In Nemo, GF2⁶'s primitive polynomial is p(z) = z⁶ + z⁴ + z³ + z + 1 while in https://web.ntpu.edu.tw/~yshan/BCH_code.pdf, it's p(z) = z⁶ + z + 1, both are irreducible primitive polynomials.
    @test generator_polynomial(BCH(6, 1)) == x ^ 6 + x ^ 4 + x ^ 3 + x + 1
    @test generator_polynomial(BCH(6, 2)) == (x ^ 6 + x ^ 4 + x ^ 3 + x + 1) * (1 + x ^ 2 + x ^ 4 + x ^ 5 + x ^ 6)
    @test generator_polynomial(BCH(6, 3)) == generator_polynomial(BCH(6, 2)) * (1 + x + x ^ 6)
    @test generator_polynomial(BCH(6, 4)) == generator_polynomial(BCH(6, 3)) * (1 + x ^ 3 + x ^ 6)
    @test generator_polynomial(BCH(6, 5)) == generator_polynomial(BCH(6, 4)) * (1 + x + x ^ 3)
    @test generator_polynomial(BCH(6, 6)) == generator_polynomial(BCH(6, 5)) * (1 + x + x ^ 2 + x ^ 5 + x ^ 6)
    @test generator_polynomial(BCH(6, 7)) == generator_polynomial(BCH(6, 6)) * (1 + x + x ^ 4 + x ^ 5 + x ^ 6)
    @test generator_polynomial(BCH(6, 10)) == generator_polynomial(BCH(6, 7)) * (1 + x + x ^ 2 + x ^ 4 + x ^ 6)
    @test generator_polynomial(BCH(6, 11)) == generator_polynomial(BCH(6, 10)) * (1 + x + x ^ 2)
    @test generator_polynomial(BCH(6, 13)) == generator_polynomial(BCH(6, 11)) * (1 + x ^ 5 + x ^ 6)  
    @test generator_polynomial(BCH(6, 15)) == generator_polynomial(BCH(6, 13)) * (1 + x ^ 2 + x ^ 3)
    
    test_cases = [(6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7), (6, 10), (6, 11), (6, 13), (6, 15)]
    results = [57 51 45 39 36 30 24 18 16 10 7]
    i = 0
    for (m, t) in test_cases
        i += 1
        @test code_k(BCH(m, t)) == results[i]
        @test check_designed_distance(parity_checks(BCH(m, t)), t) == true
    end
end
