@testitem "ECC 2BGA Table 2 via Presentation of Cyclic Groups" begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using JuMP
        using HiGHS
        using Nemo
        using Nemo: FqFieldElem
        using Hecke: group_algebra, GF, abelian_group, gens, quo, one, GroupAlgebra
        using QuantumClifford.ECC
        using QuantumClifford.ECC: code_k, code_n, two_block_group_algebra_codes, twobga_from_fp_group, 
        using Oscar: free_group, small_group_identification, describe, order, FPGroupElem, FPGroup, FPGroupElem

        @testset "Reproduce Table 2 lin2024quantum" begin
            # codes taken from Appendix C, Table 2 of [lin2024quantum](@cite)

            # [[16, 2, 4]]
            m = 4
            F = free_group(["x", "s"])
            x, s = gens(F)
            G, = quo(F, [x^m, s^2, x * s * x^-1 * s^-1])
            F2G = group_algebra(GF(2), G)
            x, s = gens(G)
            a_elts = [one(G), x]
            b_elts = [one(G), x, s, x^2, s*x, s*x^3]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 16 && code_k(c) == 2
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 4
            describe(small_group(2*m, 2))
            @test small_group_identification(G) == (order(G), 2)

            # [[16, 4, 4]]
            b_elts = [one(G), x, s, x^2, s*x, x^3]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 16 && code_k(c) == 4
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 4
            describe(small_group(2*m, 2))
            @test small_group_identification(G) == (order(G), 2)

            # [[16, 8, 2]]
            a_elts = [one(G), s]
            b_elts = [one(G), x, s, x^2, s*x, x^2]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 16 && code_k(c) == 8
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 2
            describe(small_group(2*m, 2)
            @test small_group_identification(G) == (order(G), 2)

            # [[24, 4, 5]]
            m = 6
            F = free_group(["x", "s"])
            x, s = gens(F)
            G, = quo(F, [x^m, s^2, x * s * x^-1 * s^-1])
            F2G = group_algebra(GF(2), G)
            x, s = gens(G)
            a_elts = [one(G), x]
            b_elts = [one(G), x^3, s, x^4, x^2, s*x]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 24 && code_k(c) == 4
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 5
            describe(small_group(2*m, 5))
            @test small_group_identification(G) == (order(G), 5)

            # [[24, 12, 2]]
            a_elts = [one(G), x^3]
            b_elts = [one(G), x^3, s, x^4, s * x^3, x]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 24 && code_k(c) == 12
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 2
            describe(small_group(2*m, 5))
            @test small_group_identification(G) == (order(G), 5)

            # [[32, 8, 4]]
            m = 8
            F = free_group(["x", "s"])
            x, s = gens(F)
            G, = quo(F, [x^m, s^2, x * s * x^-1 * s^-1])
            F2G = group_algebra(GF(2), G)
            x, s = gens(G)
            a_elts = [one(G), x^6]
            b_elts = [one(G), s * x^7, s * x^4, x^6, s * x^5, s * x^2]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 32 && code_k(c) == 8
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 4
            describe(small_group(2*m, 5))
            @test small_group_identification(G) == (order(G), 5)

            # [[32, 16, 2]]
            a_elts = [one(G), s * x^4]
            b_elts = [one(G), s * x^7, s * x^4, x^6, x^3, s * x^2]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 32 && code_k(c) == 16
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 2
            describe(small_group(2*m, 5))
            @test small_group_identification(G) == (order(G), 5)

            # [[40, 4, 8]]
            m = 10
            F = free_group(["x", "s"])
            x, s = gens(F)
            G, = quo(F, [x^m, s^2, x * s * x^-1 * s^-1])
            F2G = group_algebra(GF(2), G)
            x, s = gens(G)
            a_elts = [one(G), x]
            b_elts = [one(G), x^5, x^6, s * x^6, x^7, s * x^3]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 40 && code_k(c) == 4
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 8
            describe(small_group(2*m, 5))
            @test small_group_identification(G) == (order(G), 5)

            # [[40, 8, 5]]
            a_elts = [one(G), x^6]
            b_elts = [one(G), x^5, s, x^6 , x, s * x^2]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 40 && code_k(c) == 8
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 5
            describe(small_group(2*m, 5))
            @test small_group_identification(G) == (order(G), 5)

            # [[40, 20, 2]]
            a_elts = [one(G), x^5]
            b_elts = [one(G), x^5, s, x^6, s * x^5, x]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 40 && code_k(c) == 20
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 2
            describe(small_group(2*m, 5))
            @test small_group_identification(G) == (order(G), 5)

            # [[48, 8, 6]]
            m = 12
            F = free_group(["x", "s"])
            x, s = gens(F)
            G, = quo(F, [x^m, s^2, x * s * x^-1 * s^-1])
            F2G = group_algebra(GF(2), G)
            x, s = gens(G)
            a_elts = [one(G), s * x^10]
            b_elts = [one(G), x^3, s * x^6, x^4, x^7, x^8]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 48 && code_k(c) == 8
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 6
            describe(small_group(2*m, 9))
            @test small_group_identification(G) == (order(G), 9)

            # [[48, 12, 4]]
            a_elts = [one(G), x^3]
            b_elts = [one(G), x^3, s * x^6, x^4, s * x^9, x^7]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 48 && code_k(c) == 12
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 4
            describe(small_group(2*m, 9))
            @test small_group_identification(G) == (order(G), 9)

            # [[48, 16, 3]]
            a_elts = [one(G), x^4]
            b_elts = [one(G), x^3, s * x^6, x^4, x^7, s * x^10]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 48 && code_k(c) == 16
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 3
            describe(small_group(2*m, 9))
            @test small_group_identification(G) == (order(G), 9)

            # [[48, 24, 2]]
            a_elts = [one(G),  s * x^6]
            b_elts = [one(G), x^3, s * x^6, x^4, s * x^9, s * x^10]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C12 x C2"
            @test code_n(c) == 48 && code_k(c) == 24
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 2
            describe(small_group(2*m, 9))
            @test small_group_identification(G) == (order(G), 9)

            # [[56, 8, 7]]
            m = 14
            F = free_group(["x", "s"])
            x, s = gens(F)
            G, = quo(F, [x^m, s^2, x * s * x^-1 * s^-1])
            F2G = group_algebra(GF(2), G)
            x, s = gens(G)
            a_elts = [one(G), x^8]
            b_elts = [one(G), x^7, s, x^8, x^9, s * x^4]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 56 && code_k(c) == 8
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 7
            describe(small_group(2*m, 4))
            @test small_group_identification(G) == (order(G), 4)

            # [[56, 28, 2]]
            a_elts = [one(G), x^7]
            b_elts = [one(G), x^7, s, x^8, s * x^7, x]
            twobga_from_fp_group(a_elts, b_elts, F2G)
            @test order(G) == 2*m
            @test describe(G) == "C$m x C2"
            @test code_n(c) == 56 && code_k(c) == 28
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS) == 2
            describe(small_group(2*m, 4))
            @test small_group_identification(G) == (order(G), 4)
        end
    end
end
