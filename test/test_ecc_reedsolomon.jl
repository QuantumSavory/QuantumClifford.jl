using Test
using Nemo
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedSolomon

@testset "Test ReedSolomon(n, k, j) Matrix Degree" begin
    #RS (7, 3) code
    for n in 7:7
        for k in 3:3
            for j in 2:500
                deg = ceil(Int, log2(n + 1))
                t = div(n - k, 2)
                GF2ͬ, a = finite_field(2, deg, "a")
                P, x = GF2ͬ[:x]
                poly = x - a^j
                for i in 1:(2*t - 1)
                   poly *= (x - a^(j + i))
                end
                @test degree(poly) - div(deg, 2) == n - k - 1
            end
        end
    end

    # RS(255, 223) code
    for n in 255:255
        for k in 223:223
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
end