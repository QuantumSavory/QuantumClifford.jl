using Test
using LinearAlgebra
using Combinatorics
using Nemo: finite_field, GF, FqFieldElem, coeff, is_zero, degree
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedSolomon

function designed_distance(matrix, k, n, r)
    for row in eachrow(matrix)
        count = sum(row)
        # Employing 3-level quantization of the channel bits and erasing entire symbols if any of their constituent bits are erased can improve the performance of RS codes. Refer to chapter: 07, section: 03, pages: 172 to 175 [tomlinson2017error](@cite).
        if count >= 2 ^ (r + 1) - 3 - k
            return true
        end
    end
    return false
end

@testset "Testing Shortened and Expanded Maximum Distance Separable (MDS) Reed Solomon codes's binary parity check matrices" begin
    n_cases = [31, 63, 127, 255]
    for n in n_cases
        r = ceil(Int, log2(n + 1))
        for k in 5:15
        @test rank(parity_checks(ReedSolomon(n, k))) == r * k
        @test designed_distance(parity_checks(ReedSolomon(n, k)), k, n, r) == true
        end
    end
end

@testset "Test ReedSolomon(n, k) generator polynomial g(x) degree" begin
    # RS(7, 3), RS(15, 9), RS(255, 223), RS(160, 128), RS(255, 251), (255, 239) and (255, 249) codes
    test_cases = [(7, 3), (15, 9), (225, 223), (160, 128), (255, 251), (255, 239), (255, 249)]
    for (n, k) in test_cases
        deg = ceil(Int, log2(n + 1))
        t = div(n - k, 2)
        GF2ͬ, a = finite_field(2, deg, "a")
        P, x = GF2ͬ[:x]
        pzeros = 2*t
        poly = x - a^pzeros
        for i in 1:(2*t - 1)
           poly *= (x - a^(pzeros + i))
        end
        @test degree(poly) == n - k
    end
end
