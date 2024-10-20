@testitem "ECC 2BGA abelian and non-abelian groups via group presentation" begin
    using Nemo: FqFieldElem
    using Hecke: group_algebra, GF, abelian_group, gens, quo, one, GroupAlgebra
    using QuantumClifford.ECC
    using QuantumClifford.ECC: code_k, code_n, two_block_group_algebra_codes
    using Oscar: free_group, small_group_identification, describe, order, FPGroupElem, FPGroup, FPGroupElem

    function get_code(a_elts::Vector{FPGroupElem}, b_elts::Vector{FPGroupElem}, F2G::GroupAlgebra{FqFieldElem, FPGroup, FPGroupElem})
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a,b)
        return c
    end

    @testset "Reproduce Table 1 Block 1" begin
        # [[72, 8, 9]]
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^36])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^28]
        b_elts = [one(G), r, r^18, r^12, r^29, r^14]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C36"
        @test small_group_identification(G) == (36, 2)
        @test code_n(c) == 72 && code_k(c) == 8
    end

    @testset "Reproduce Table 1 Block 2" begin
        # [[72, 8, 9]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^4, r^9, s^(-1)*r*s*r])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r]
        b_elts = [one(G), s, r^6, s^3 * r, s * r^7, s^3 * r^5]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C9 : C4"
        @test small_group_identification(G) == (36, 1)
        @test code_n(c) == 72 && code_k(c) == 8

        # [[80, 8, 10]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^5, r^8, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), s*r^4]
        b_elts = [one(G), r, r^2, s, s^3 * r, s^2 * r^6]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C5 : C8"
        @test small_group_identification(G) == (40, 1)
        @test code_n(c) == 80 && code_k(c) == 8
    end

    @testset "Reproduce Table 1 Block 3" begin
        # [[54, 6, 9]]
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^27])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r, r^3, r^7]
        b_elts = [one(G), r, r^12, r^19]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C27"
        @test small_group_identification(G) == (27, 1)
        @test code_n(c) == 54 && code_k(c) == 6

        # [[60, 6, 10]]
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^30])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^10, r^6, r^13]
        b_elts = [one(G), r^25, r^16, r^12]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C30"
        @test small_group_identification(G) == (30, 4)
        @test code_n(c) == 60 && code_k(c) == 6

        # [[70, 8, 10]]
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^35])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^15, r^16, r^18]
        b_elts = [one(G), r, r^24, r^27]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C35"
        @test small_group_identification(G) == (35, 1)
        @test code_n(c) == 70 && code_k(c) == 8

        # [[72, 8, 10]]
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^36])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^9, r^28, r^31]
        b_elts = [one(G), r, r^21, r^34]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C36"
        @test small_group_identification(G) == (36, 2)
        @test code_n(c) == 72 && code_k(c) == 8

        # [[72, 10, 9]]
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^36])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^9, r^28, r^13]
        b_elts = [one(G), r, r^3, r^22]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C36"
        @test small_group_identification(G) == (36, 2)
        @test code_n(c) == 72 && code_k(c) == 10
    end

    @testset "Reproduce Table 1 Block 4" begin
        # [[72, 8, 9]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^4, r^9, s^(-1)*r*s*r])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), s, r, s*r^6]
        b_elts = [one(G), s^2*r, s^2*r^6, r^2]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C9 : C4"
        @test small_group_identification(G) == (36, 1)
        @test code_n(c) == 72 && code_k(c) == 8

        # [[80, 8, 10]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^5, r^8, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r, s, s^3*r^5]
        b_elts = [one(G), r^2, s*r^4, s^3*r^2]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C5 : C8"
        @test small_group_identification(G) == (40, 1)
        @test code_n(c) == 80 && code_k(c) == 8

        # [[96, 8, 12]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^3, r^16, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r, s, r^14]
        b_elts = [one(G), r^2, s*r^4, r^11]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C3 : C16"
        @test small_group_identification(G) == (48, 1)
        @test code_n(c) == 96 && code_k(c) == 6

        # [[84, 10, 9]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^3, r^14, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^7, r^8, s*r^10]
        b_elts = [one(G), s, r^5, s^2*r^13]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C7 x S3"
        @test small_group_identification(G) == (42, 3)
        @test code_n(c) == 84 && code_k(c) == 10

        # [[96, 6, 12]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^4, r^12, s^(-1)*r*s*r])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), s, r^9, s * r]
        b_elts = [one(G), s^2 * s^9, r^7, r^2]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C12 : C4"
        @test small_group_identification(G) == (48, 13)
        @test code_n(c) == 96 && code_k(c) == 6

        # [[96, 12, 10]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^6, r^8, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r, s^3 * r^2, s^2 * r^3]
        b_elts = [one(G), r, s^4 * r^6, s^5 * r^3]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C2 x (C3 : C8)"
        @test small_group_identification(G) == (48, 9)
        @test code_n(c) == 96 && code_k(c) == 12
    end
end
