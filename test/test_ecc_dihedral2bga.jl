@testitem "ECC 2BGA Reproduce Table 3 of arXiv:2306.16400" tags=[:ecc, :ecc_bespoke_checks, :oscar_required] begin
    using Nemo
    using Nemo: FqFieldElem
    using Hecke: group_algebra, GF, abelian_group, gens, quo, one, GroupAlgebra
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: code_k, code_n, two_block_group_algebra_codes, DistanceMIPAlgorithm
    using Oscar
    using Oscar: free_group, small_group_identification, describe, order, FPGroupElem, FPGroup, FPGroupElem
    using HiGHS
    using JuMP

    @testset "Reproduce Table 3 of arXiv:2306.16400" begin
        # [[24, 8, 3]]
        m = 6
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^4]
        b_elts = [one(G), s*r^4, r^3, r^4, s*r^2, r]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 24 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test small_group_identification(G) == (order(G), 4)

        # [[24, 12, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^3]
        b_elts = [one(G), s*r, r^3, r^4, s*r^4, r]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 24 && code_k(c) == 12
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 2
        @test small_group_identification(G) == (order(G), 4)

        # [[32, 8, 4]]
        m = 8
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^2]
        b_elts = [one(G), s*r^5, s*r^4, r^2, s*r^7, s*r^6]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 32 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        @test small_group_identification(G) == (order(G), 7)

        # [[32, 16, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^4]
        b_elts = [one(G), s*r^3, s*r^6, r^4, s*r^7, s*r^2]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 32 && code_k(c) == 16
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 2
        @test small_group_identification(G) == (order(G), 7)

        # [[36, 12, 3]]
        m = 9
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^3]
        b_elts = [one(G), s, r, r^3, s*r^3, r^4]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 36 && code_k(c) == 12
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test small_group_identification(G) == (order(G), 1)

        # [[40, 8, 5]]
        m = 10
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^2]
        b_elts = [one(G), s*r^4, r^5, r^2, s*r^6, r]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 40 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5
        @test small_group_identification(G) == (order(G), 4)

        # [[40, 20, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^5]
        b_elts = [one(G), s*r^2, r^5, r^6, s*r^7, r]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 40 && code_k(c) == 20
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 2
        @test small_group_identification(G) == (order(G), 4)

        # [[48, 8, 6]]
        m = 12
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^10]
        b_elts = [one(G), s*r^8, r^9, r^4, s*r^2, r^5]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 48 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6
        @test small_group_identification(G) == (order(G), 6)

        # [[48, 12, 4]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^3]
        b_elts = [one(G), s*r^7, r^3, r^4, s*r^10, r^7]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 48 && code_k(c) == 12
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        @test small_group_identification(G) == (order(G), 6)

        # [[48, 16, 3]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^8]
        b_elts = [one(G), s*r^8, r^9, r^8, s*r^4, r^5]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 48 && code_k(c) == 16
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test small_group_identification(G) == (order(G), 6)

        # [[48, 24, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^6]
        b_elts = [one(G), s*r^11, r^6, s*r^5, r, r^7]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 48 && code_k(c) == 24
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 2
        @test small_group_identification(G) == (order(G), 6)

        # [[56, 8, 7]]
        m = 14
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^4]
        b_elts = [one(G), s*r^11, r^7, s*r^5, r^12, r^9]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 56 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 7
        @test small_group_identification(G) == (order(G), 3)

        # [[56, 28, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^7]
        b_elts = [one(G), s*r^2, r^7, r^8, s*r^9, r]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 56 && code_k(c) == 28
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 2
        @test small_group_identification(G) == (order(G), 3)

        # [[60, 12, 5]]
        m = 15
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^12]
        b_elts = [one(G), s*r^14, r^5, r^12, s*r^11, r^14]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 60 && code_k(c) == 12
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5
        @test small_group_identification(G) == (order(G), 3)

        # [[60, 20, 3]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^5]
        b_elts = [one(G), s*r^13, r^5, r^12, s*r^3, r^2]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 60 && code_k(c) == 20
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test small_group_identification(G) == (order(G), 3)

        # [[64, 8, 8]]
        m = 16
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^6]
        b_elts = [one(G), s*r^12, s*r^9, r^6, s, s*r]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 64 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8
        @test small_group_identification(G) == (order(G), 18)

        # [[64, 16, 4]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^4]
        b_elts = [one(G), s*r^10, s*r^3, r^4, s*r^14, s*r^7]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 64 && code_k(c) == 16
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        @test small_group_identification(G) == (order(G), 18)

        # [[64, 32, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^m, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^8]
        b_elts = [one(G), s*r^11, s*r^12, r^8, s*r^3, s*r^4]
        c = twobga_from_fp_group(a_elts, b_elts, F2G)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 64 && code_k(c) == 32
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 2
        @test small_group_identification(G) == (order(G), 18)
    end
end
