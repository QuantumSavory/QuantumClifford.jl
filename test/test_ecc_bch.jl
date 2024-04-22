using Test
using Nemo: ZZ, residue_ring, matrix, finite_field, GF, minpoly, coeff, lcm, FqPolyRingElem, FqFieldElem, is_zero
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, BCH

@testset "Test BCH(n, t) generator polynomial g(x) universal property: mod(x^n - 1, g(x)) == 0" begin
    n_cases = [7, 15, 31, 63, 127, 255]
    for n in n_cases
        #Testing all 1 Bit, 2 Bit, 3 Bit and 4 Bit BCH codes for n_cases
        for t in 1:4
            gx = generator_polynomial(BCH(n, t))
            GF2x, x = GF(2)["x"] 
            @test mod(x^n - 1, gx) == 0
        end
    end
end

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

@testset "Test 1 Bit and 2 Bit BCH codes designed distance from binary parity check matrix H" begin
    n_cases = [7, 15, 31, 63, 127, 255]
    for n in n_cases
        for t in 1:2
            r = ceil(Int, log2(n + 1))
            H = parity_checks(BCH(n, t))
            @test designed_distance(H) >= 2 * t + 1 || designed_distance(H) == 2 * t
        end
    end
end
