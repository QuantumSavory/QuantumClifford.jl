using Test
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, BCH, generator_polynomial
using Nemo: ZZ, residue_ring, matrix, finite_field, GF, minpoly, coeff, lcm, FqPolyRingElem, FqFieldElem, is_zero, degree, defining_polynomial, is_irreducible 
include("utils_test_ecc.jl")

"""
- To prove that `t`-bit error correcting BCH code indeed has minimum distance of at least `2 * t + 1`, it is shown that no `2 * t` or fewer columns of its binary parity check matrix `H` sum to zero. A formal mathematical proof can be found on pages 168 and 169 of Ch6 of Error Control Coding by Lin, Shu and Costello, Daniel. 
- The parameter `2 * t + 1` is usually called the designed distance of the `t`-bit error correcting BCH code.
"""

@testset "Testing designed distance of BCH codes" begin
    m_cases = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    for m in m_cases  
        for t in rand(1:m, 2)
            d = 2 * t
            @test check_designed_distance(parity_checks(BCH(m, t)), m, t, d, 0, 0) == true
        end
    end
end

@testset "Testing properties of BCH codes" begin
    m_cases = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    for m in m_cases
        n = 2 ^ m - 1
        for t in rand(1:m, 2)
            H = parity_checks(BCH(m, t))
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

    #example taken from Ch6 of Error Control Coding by Lin, Shu and Costello, Daniel
    @test parity_checks(BCH(4, 2))  ==    [1  0  0  0  1  0  0  1  1  0  1  0  1  1  1;
                                           0  1  0  0  1  1  0  1  0  1  1  1  1  0  0;
                                           0  0  1  0  0  1  1  0  1  0  1  1  1  1  0;
                                           0  0  0  1  0  0  1  1  0  1  0  1  1  1  1;
                                           1  0  0  0  1  1  0  0  0  1  1  0  0  0  1;
                                           0  0  0  1  1  0  0  0  1  1  0  0  0  1  1;
                                           0  0  1  0  1  0  0  1  0  1  0  0  1  0  1;
                                           0  1  1  1  1  0  1  1  1  1  0  1  1  1  1]
    
    # Examples taken from https://web.ntpu.edu.tw/~yshan/BCH_code.pdf.
    GF2x, x = GF(2)["x"]
    GF2⁴, a = finite_field(2, 4, "a")
    GF2⁶, b = finite_field(2, 6, "b")
    @test defining_polynomial(GF2x, GF2⁴) ==  x ^ 4 + x + 1
    @test is_irreducible(defining_polynomial(GF2x, GF2⁴)) == true
    @test generator_polynomial(BCH(4, 2)) == x ^ 8 + x ^ 7 +  x ^ 6 +  x ^ 4 + 1
    @test generator_polynomial(BCH(4, 3)) == x ^ 10 + x ^ 8 +  x ^ 5 +  x ^ 4 +  x ^ 2 + x + 1
 
    # Nemo.jl uses [Conway polynomial](https://en.wikipedia.org/wiki/Conway_polynomial_(finite_fields)), a standard way to represent the primitive polynomial for finite Galois fields `GF(pᵐ)` of degree `m`, where `p` is a prime number. 
    # The `GF(2⁶)`'s Conway polynomial is `p(z) = z⁶ + z⁴ + z³ + z + 1`. In contrast, the polynomial given in https://web.ntpu.edu.tw/~yshan/BCH_code.pdf is `p(z) = z⁶ + z + 1`. Because both polynomials are irreducible, they are also primitive polynomials for `GF(2⁶)`.
    
    test_cases = [(6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7), (6, 10), (6, 11), (6, 13), (6, 15)]
    @test defining_polynomial(GF2x, GF2⁶) == x ^ 6 + x ^ 4 + x ^ 3 + x + 1
    @test is_irreducible(defining_polynomial(GF2x, GF2⁶)) == true    
    for i in 1:length(test_cases)
        m, t = test_cases[i]
        if t == 1
            @test generator_polynomial(BCH(m, t)) == defining_polynomial(GF2x, GF2⁶)
        else
            prev_t = test_cases[i - 1][2]
            @test generator_polynomial(BCH(m, t)) == generator_polynomial(BCH(m, prev_t)) * minpoly(GF2x, b ^ (t + prev_t - (t - prev_t - 1)))
        end
    end
    
    results = [57 51 45 39 36 30 24 18 16 10 7]
    for (result, (m, t)) in zip(results, test_cases)
        d = 2 * t
        @test code_k(BCH(m, t)) == result
        @test check_designed_distance(parity_checks(BCH(m, t)), m, t, d, 0, 0) == true
    end
end
