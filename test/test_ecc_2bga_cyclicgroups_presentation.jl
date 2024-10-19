@testitem "ECC 2BGA Table 2 via Presentation of Cyclic Groups" begin
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

    @testset "Reproduce Table 2 lin2024quantum" begin
        # codes taken from Appendix C, Table 2 of [lin2024quantum](@cite)

        # [[16, 2, 4]]
        # m = 4
        F = free_group(["x", "s"])
        x, s = gens(F)
        G, = quo(F, [x^4, s^2, x * s * x^-1 * s^-1])
        F2G = group_algebra(GF(2), G)
        x, s = gens(G)
        a_elts = [one(G), x]
        b_elts = [one(G), x, s, x^2, s*x, s*x^3]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C4 x C2"
        @test order(G) == 8
        @test small_group_identification(G) == (8, 2)
        @test code_n(c) == 16 && code_k(c) == 2

        # [[16, 4, 4]]
        b_elts = [one(G), x, s, x^2, s*x, x^3]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C4 x C2"
        @test order(G) == 8
        @test small_group_identification(G) == (8, 2)
        @test code_n(c) == 16 && code_k(c) == 4

        # [[16, 8, 2]]
        a_elts = [one(G), s]
        b_elts = [one(G), x, s, x^2, s*x, x^2]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C4 x C2"
        @test order(G) == 8
        @test small_group_identification(G) == (8, 2)
        @test code_n(c) == 16 && code_k(c) == 8

        # [[24, 4, 5]]
        # m = 6
        F = free_group(["x", "s"])
        x, s = gens(F)
        G, = quo(F, [x^6, s^2, x * s * x^-1 * s^-1])
        F2G = group_algebra(GF(2), G)
        x, s = gens(G)
        a_elts = [one(G), x]
        b_elts = [one(G), x^3, s, x^4, x^2, s*x]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C6 x C2"
        @test order(G) == 12
        @test small_group_identification(G) == (12, 5)
        @test code_n(c) == 24 && code_k(c) == 4

        # [[24, 12, 2]]
        a_elts = [one(G), x^3]
        b_elts = [one(G), x^3, s, x^4, s * x^3, x]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C6 x C2"
        @test order(G) == 12
        @test small_group_identification(G) == (12, 5)
        @test code_n(c) == 24 && code_k(c) == 12

        # [[32, 8, 4]]
        # m = 8
        F = free_group(["x", "s"])
        x, s = gens(F)
        G, = quo(F, [x^8, s^2, x * s * x^-1 * s^-1])
        F2G = group_algebra(GF(2), G)
        x, s = gens(G)
        a_elts = [one(G), x^6]
        b_elts = [one(G), s * x^7, s * x^4, x^6, s * x^5, s * x^2]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C8 x C2"
        @test order(G) == 16
        @test small_group_identification(G) == (16, 5)
        @test code_n(c) == 32 && code_k(c) == 8

        # [[32, 16, 2]]
        a_elts = [one(G), s * x^4]
        b_elts = [one(G), s * x^7, s * x^4, x^6, x^3, s * x^2]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C8 x C2"
        @test order(G) == 16
        @test small_group_identification(G) == (16, 5)
        @test code_n(c) == 32 && code_k(c) == 16

        # [[40, 4, 8]]
        # m = 10
        F = free_group(["x", "s"])
        x, s = gens(F)
        G, = quo(F, [x^10, s^2, x * s * x^-1 * s^-1])
        F2G = group_algebra(GF(2), G)
        x, s = gens(G)
        a_elts = [one(G), x]
        b_elts = [one(G), x^5, x^6, s * x^6, x^7, s * x^3]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C10 x C2"
        @test order(G) == 20
        @test small_group_identification(G) == (20, 5)
        @test code_n(c) == 40 && code_k(c) == 4

        # [[40, 8, 5]]
        a_elts = [one(G), x^6]
        b_elts = [one(G), x^5, s, x^6 , x, s * x^2]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C10 x C2"
        @test order(G) == 20
        @test small_group_identification(G) == (20, 5)
        @test code_n(c) == 40 && code_k(c) == 8

        # [[40, 20, 2]]
        a_elts = [one(G), x^5]
        b_elts = [one(G), x^5, s, x^6, s * x^5, x]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C10 x C2"
        @test order(G) == 20
        @test small_group_identification(G) == (20, 5)
        @test code_n(c) == 40 && code_k(c) == 20

        # [[48, 8, 6]]
        # m = 12
        F = free_group(["x", "s"])
        x, s = gens(F)
        G, = quo(F, [x^12, s^2, x * s * x^-1 * s^-1])
        F2G = group_algebra(GF(2), G)
        x, s = gens(G)
        a_elts = [one(G), s * x^10]
        b_elts = [one(G), x^3, s * x^6, x^4, x^7, x^8]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C12 x C2"
        @test order(G) == 24
        @test small_group_identification(G) == (24, 9)
        @test code_n(c) == 48 && code_k(c) == 8

        # [[48, 12, 4]]
        a_elts = [one(G), x^3]
        b_elts = [one(G), x^3, s * x^6, x^4, s * x^9, x^7]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C12 x C2"
        @test order(G) == 24
        @test small_group_identification(G) == (24, 9)
        @test code_n(c) == 48 && code_k(c) == 12

        # [[48, 16, 3]]
        a_elts = [one(G), x^4]
        b_elts = [one(G), x^3, s * x^6, x^4, x^7, s * x^10]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C12 x C2"
        @test order(G) == 24
        @test small_group_identification(G) == (24, 9)
        @test code_n(c) == 48 && code_k(c) == 16

        # [[48, 24, 2]]
        a_elts = [one(G),  s * x^6]
        b_elts = [one(G), x^3, s * x^6, x^4, s * x^9, s * x^10]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C12 x C2"
        @test order(G) == 24
        @test small_group_identification(G) == (24, 9)
        @test code_n(c) == 48 && code_k(c) == 24

        # [[56, 8, 7]]
        # m = 14
        F = free_group(["x", "s"])
        x, s = gens(F)
        G, = quo(F, [x^14, s^2, x * s * x^-1 * s^-1])
        F2G = group_algebra(GF(2), G)
        x, s = gens(G)
        a_elts = [one(G), x^8]
        b_elts = [one(G), x^7, s, x^8, x^9, s * x^4]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C14 x C2"
        @test order(G) == 28
        @test small_group_identification(G) == (28, 4)
        @test code_n(c) == 56 && code_k(c) == 8

        # [[56, 28, 2]]
        a_elts = [one(G), x^7]
        b_elts = [one(G), x^7, s, x^8, s * x^7, x]
        c = get_code(a_elts, b_elts, F2G)
        @test describe(G) == "C14 x C2"
        @test order(G) == 28
        @test small_group_identification(G) == (28, 4)
        @test code_n(c) == 56 && code_k(c) == 28
    end
end
