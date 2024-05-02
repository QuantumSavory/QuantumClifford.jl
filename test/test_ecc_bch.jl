using Test
using Nemo: ZZ, residue_ring, matrix, finite_field, GF, minpoly, coeff, lcm, FqPolyRingElem, FqFieldElem, is_zero, degree, matrix
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, BCH

function designed_distance(matrix)
    mindist = Inf 
    for row in eachrow(matrix)
        count = sum(row)
        if count < mindist
            mindist = count
        end
    end
    return mindist
end

@testset "Testing properties of BCH codes for all its instances" begin
    m_cases = [3, 4, 5, 6, 7, 8, 9, 10]
    for m in m_cases
        n = 2 ^ m - 1
        for t in 1:2
            H = parity_checks(BCH(m, t))
            @test designed_distance(H) >= 2 * t
            # n - k == degree of generator polynomial, gx
            mat = matrix(GF(2), parity_checks(BCH(m, t)))
            computed_rank = rank(mat)
            @test computed_rank == degree(generator_polynomial(BCH(m, t)))
        end
        # BCH code is cyclic as its generator polynomial, gx divides xⁿ - 1, so mod (xⁿ - 1, gx) = 0 
        for t in 1:2
            gx = generator_polynomial(BCH(m, t))
            GF2x, x = GF(2)["x"] 
            @test mod(x^n - 1, gx) == 0
        end
    end

    #example taken from [https://web.ntpu.edu.tw/~yshan/BCH_code.pdf]    
    @test parity_checks(BCH(4, 2))  ==    [1  0  0  0  1  0  0  1  1  0  1  0  1  1  1;
                                           0  1  0  0  1  1  0  1  0  1  1  1  1  0  0;
                                           0  0  1  0  0  1  1  0  1  0  1  1  1  1  0;
                                           0  0  0  1  0  0  1  1  0  1  0  1  1  1  1;
                                           1  0  0  0  1  1  0  0  0  1  1  0  0  0  1;
                                           0  0  0  1  1  0  0  0  1  1  0  0  0  1  1;
                                           0  0  1  0  1  0  0  1  0  1  0  0  1  0  1;
                                           0  1  1  1  1  0  1  1  1  1  0  1  1  1  1]

    GF2ͬ, a = finite_field(2, 4, "a")
    GF2x, x = GF(2)["x"]
    #examples taken from [https://web.ntpu.edu.tw/~yshan/BCH_code.pdf]   
    @test generator_polynomial(BCH(4, 2)) == x ^ 8 + x ^ 7 +  x ^ 6 +  x ^ 4 + 1
    @test generator_polynomial(BCH(4, 3)) == x ^ 10 + x ^ 8 +  x ^ 5 +  x ^ 4 +  x ^ 2 + x + 1

end