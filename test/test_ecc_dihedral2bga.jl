@testitem "ECC 2BGA lin2024quantum" begin
    using Hecke
    using Hecke: group_algebra, GF, abelian_group, gens, quo, one
    using QuantumClifford.ECC: LPCode, code_k, code_n
    using Oscar: free_group, small_group_identification, describe, order

    @testset "Reproduce Table 3 of lin2024quantum" begin
        # [[24, 8, 3]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^6, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^4]
        b_elts = [one(G), s*r^4, r^3, r^4, s*r^2, r]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D12"
        @test order(G) == 12
        @test small_group_identification(G) == (12, 4)
        @test code_n(c) == 24 && code_k(c) == 8

        # [[24, 12, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^6, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^3]
        b_elts = [one(G), s*r, r^3, r^4, s*r^4, r]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D12"
        @test order(G) == 12
        @test small_group_identification(G) == (12, 4)
        @test code_n(c) == 24 && code_k(c) == 12

        # [[32, 8, 4]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^8, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^2]
        b_elts = [one(G), s*r^5, s*r^4, r^2, s*r^7, s*r^6]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D16"
        @test order(G) == 16
        @test small_group_identification(G) == (16, 7)
        @test code_n(c) == 32 && code_k(c) == 8

        # [[32, 16, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^8, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^4]
        b_elts = [one(G), s*r^3, s*r^6, r^4, s*r^7, s*r^2]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D16"
        @test order(G) == 16
        @test small_group_identification(G) == (16, 7)
        @test code_n(c) == 32 && code_k(c) == 16

        # [[36, 12, 3]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^9, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^3]
        b_elts = [one(G), s, r, r^3, s*r^3, r^4]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D18"
        @test order(G) == 18
        @test small_group_identification(G) == (18, 1)
        @test code_n(c) == 36 && code_k(c) == 12

        # [[40, 8, 5]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^10, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^2]
        b_elts = [one(G), s*r^4, r^5, r^2, s*r^6, r]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D20"
        @test order(G) == 20
        @test small_group_identification(G) == (20, 4)
        @test code_n(c) == 40 && code_k(c) == 8

        # [[40, 20, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^10, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^5]
        b_elts = [one(G), s*r^2, r^5, r^6, s*r^7, r]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D20"
        @test order(G) == 20
        @test small_group_identification(G) == (20, 4)
        @test code_n(c) == 40 && code_k(c) == 20

        # [[48, 8, 6]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^12, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^10]
        b_elts = [one(G), s*r^8, r^9, r^4, s*r^2, r^5]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D24"
        @test order(G) == 24
        @test small_group_identification(G) == (24, 6)
        @test code_n(c) == 48 && code_k(c) == 8

        # [[48, 12, 4]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^12, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^3]
        b_elts = [one(G), s*r^7, r^3, r^4, s*r^10, r^7]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D24"
        @test order(G) == 24
        @test small_group_identification(G) == (24, 6)
        @test code_n(c) == 48 && code_k(c) == 12

        # [[48, 16, 3]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^12, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^8]
        b_elts = [one(G), s*r^8, r^9, r^8, s*r^4, r^5]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D24"
        @test order(G) == 24
        @test small_group_identification(G) == (24, 6)
        @test code_n(c) == 48 && code_k(c) == 16

        # [[48, 24, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^12, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^6]
        b_elts = [one(G), s*r^11, r^6, s*r^5, r, r^7]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D24"
        @test order(G) == 24
        @test small_group_identification(G) == (24, 6)
        @test code_n(c) == 48 && code_k(c) == 24

        # [[56, 8, 7]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^14, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^4]
        b_elts = [one(G), s*r^11, r^7, s*r^5, r^12, r^9]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D28"
        @test order(G) == 28
        @test small_group_identification(G) == (28, 3)
        @test code_n(c) == 56 && code_k(c) == 8

        # [[56, 28, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^14, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^7]
        b_elts = [one(G), s*r^2, r^7, r^8, s*r^9, r]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D28"
        @test order(G) == 28
        @test small_group_identification(G) == (28, 3)
        @test code_n(c) == 56 && code_k(c) == 28


        # [[60, 12, 5]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^15, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^12]
        b_elts = [one(G), s*r^14, r^5, r^12, s*r^11, r^14]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D30"
        @test order(G) == 30
        @test small_group_identification(G) == (30, 3)
        @test code_n(c) == 60 && code_k(c) == 12

        # [[60, 20, 3]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^15, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^5]
        b_elts = [one(G), s*r^13, r^5, r^12, s*r^3, r^2]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D30"
        @test order(G) == 30
        @test small_group_identification(G) == (30, 3)
        @test code_n(c) == 60 && code_k(c) == 20

        # [[64, 8, 8]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^16, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^6]
        b_elts = [one(G), s*r^12, s*r^9, r^6, s, s*r]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D32"
        @test order(G) == 32
        @test small_group_identification(G) == (32, 18)
        @test code_n(c) == 64 && code_k(c) == 8

        # [[64, 16, 8]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^16, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^4]
        b_elts = [one(G), s*r^10, s*r^3, r^4, s*r^14, s*r^7]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D32"
        @test order(G) == 32
        @test small_group_identification(G) == (32, 18)
        @test code_n(c) == 64 && code_k(c) == 16

        # [[64, 32, 2]]
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [r^16, s^2, (r*s)^2])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^8]
        b_elts = [one(G), s*r^11, s*r^12, r^8, s*r^3, s*r^4]
        a = sum(F2G(x) for x in a_elts)
        b = sum(F2G(x) for x in b_elts)
        c = two_block_group_algebra_codes(a, b)
        @test describe(G) == "D32"
        @test order(G) == 32
        @test small_group_identification(G) == (32, 18)
        @test code_n(c) == 64 && code_k(c) == 32
    end
end
