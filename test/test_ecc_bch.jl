using Test
using Nemo
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC

@testset "Test BCH(n, t) generator polynomial g(x) universal property: mod(x^n - 1, g(x)) == 0" begin
    n_cases = [7, 15, 31, 63, 127, 255]
    for n in n_cases
        #Testing all 1 Bit, 2 Bit, 3 Bit and 4 Bit BCH codes for n_cases
        for t in 1:4
            r = ceil(Int, log2(n + 1))
            GF2ͬ, a = finite_field(2, r, "a")
            GF2x, x = GF(2)["x"] 
            mx = FqPolyRingElem[]
            for i in 1:(2*t - 1)
                if i % 2 != 0
                    mx = [mx; minpoly(GF2x, a^i)]
                end 
            end
            gx = lcm(mx)
            @test mod(x^n - 1, gx) == 0
        end
    end
end
