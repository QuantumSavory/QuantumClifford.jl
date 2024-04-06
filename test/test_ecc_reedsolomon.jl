using Test
using Nemo
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedSolomon

@testset "Test ReedSolomon(n, k, j)  generator polynomial g(x) degree" begin
    # RS(7, 3), RS((15, 9), RS(255, 223), RS(160, 128), RS(255, 251), (255, 239) and (255, 249) codes
    test_cases = [(7, 3), (15, 9), (225, 223), (160, 128), (255, 251), (255, 239), (255, 249)]
    for (n, k) in test_cases
        for j in 2:500
            deg = ceil(Int, log2(n + 1))
            t = div(n - k, 2)
            GF2ͬ, a = finite_field(2, deg, "a")
            P, x = GF2ͬ[:x]
            poly = x - a^j
            for i in 1:(2*t - 1)
               poly *= (x - a^(j + i))
            end
            @test degree(poly) == n - k
        end 
    end
end
