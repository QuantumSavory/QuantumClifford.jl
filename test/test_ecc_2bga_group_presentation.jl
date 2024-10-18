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
        l = 36
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^l])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^28]
        b_elts = [one(G), r, r^18, r^12, r^29, r^14]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l
        @test describe(G) == "C$l"
        @test code_n(c) == 72 && code_k(c) == 8
        # Oscar.describe(Oscar.small_group(l, 2)) is C₃₆, cross-check it with G
        @test small_group_identification(G) == (order(G), 2)
    end

    @testset "Reproduce Table 1 Block 2" begin
        # [[72, 8, 9]]
        l = 9
        m = 4
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^m, r^l, s^(-1)*r*s*r])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r]
        b_elts = [one(G), s, r^6, s^3 * r, s * r^7, s^3 * r^5]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l*m
        @test describe(G) == "C$l : C$m"
        @test code_n(c) == 72 && code_k(c) == 8
        # Oscar.describe(Oscar.small_group(l*m, 1)) is C₉ x C₄, cross-check it with G
        @test small_group_identification(G) == (order(G), 1)

        # [[80, 8, 10]]
        l = 5
        m = 8
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^l, r^m, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), s*r^4]
        b_elts = [one(G), r, r^2, s, s^3 * r, s^2 * r^6]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l*m
        @test describe(G) == "C$l : C$m"
        @test code_n(c) == 80 && code_k(c) == 8
        # Oscar.describe(Oscar.small_group(l*m, 1)) is C₅ x C₈, cross-check it with G
        @test small_group_identification(G) == (order(G), 1)
    end

    @testset "Reproduce Table 1 Block 3" begin
        # [[54, 6, 9]]
        l = 27
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^l])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r, r^3, r^7]
        b_elts = [one(G), r, r^12, r^19]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l
        @test describe(G) == "C$l"
        @test code_n(c) == 54 && code_k(c) == 6
        # Oscar.describe(Oscar.small_group(l, 1)) is C₂₇, cross-check it with G
        @test small_group_identification(G) == (order(G), 1)

        # [[60, 6, 10]]
        l = 30
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^l])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^10, r^6, r^13]
        b_elts = [one(G), r^25, r^16, r^12]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l
        @test describe(G) == "C$l"
        @test code_n(c) == 60 && code_k(c) == 6
        # Oscar.describe(Oscar.small_group(l, 4)) is C₃₀, cross-check it with G
        @test small_group_identification(G) == (order(G), 4)

        # [[70, 8, 10]]
        l = 35
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^l])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^15, r^16, r^18]
        b_elts = [one(G), r, r^24, r^27]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l
        @test describe(G) == "C$l"
        @test code_n(c) == 70 && code_k(c) == 8
        # Oscar.describe(Oscar.small_group(l, 1)) is C₃₅, cross-check it with G
        @test small_group_identification(G) == (order(G), 1)

        # [[72, 8, 10]]
        l = 36
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^l])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^9, r^28, r^31]
        b_elts = [one(G), r, r^21, r^34]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l
        @test describe(G) == "C$l"
        @test code_n(c) == 72 && code_k(c) == 8
        # Oscar.describe(Oscar.small_group(l, 2)) is C₃₆, cross-check it with G
        @test small_group_identification(G) == (order(G), 2)

        # [[72, 10, 9]]
        F = free_group(["r"])
        r = gens(F)[1]
        G, = quo(F, [r^l])
        F2G = group_algebra(GF(2), G)
        r = gens(G)[1]
        a_elts = [one(G), r^9, r^28, r^13]
        b_elts = [one(G), r, r^3, r^22]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l
        @test describe(G) == "C$l"
        @test code_n(c) == 72 && code_k(c) == 10
        # Oscar.describe(Oscar.small_group(l, 2)) is C₃₆, cross-check it with G
        @test small_group_identification(G) == (order(G), 2)
    end

    @testset "Reproduce Table 1 Block 4" begin
        # [[72, 8, 9]]
        l = 9
        m = 4
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^m, r^l, s^(-1)*r*s*r])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), s, r, s*r^6]
        b_elts = [one(G), s^2*r, s^2*r^6, r^2]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l*m
        @test describe(G) == "C$l : C$m"
        @test code_n(c) == 72 && code_k(c) == 8
        # Oscar.describe(Oscar.small_group(l*m, 1)) is C₉ x C₄, cross-check it with G
        @test small_group_identification(G) == (order(G), 1)

        # [[80, 8, 10]]
        l = 5
        m = 8
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^l, r^m, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r, s, s^3*r^5]
        b_elts = [one(G), r^2, s*r^4, s^3*r^2]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l*m
        @test describe(G) == "C$l : C$m"
        @test code_n(c) == 80 && code_k(c) == 8
        # Oscar.describe(Oscar.small_group(l*m, 1)) is C₅ x C₈, cross-check it with G
        @test small_group_identification(G) == (order(G), 1)

        # [[96, 8, 12]]
        l = 3
        m = 16
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^l, r^m, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r, s, r^14]
        b_elts = [one(G), r^2, s*r^4, r^11]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l*m
        @test describe(G) == "C$l : C$m"
        @test code_n(c) == 96 && code_k(c) == 6
        # Oscar.describe(Oscar.small_group(l*m, 1)) is C₃ x C₁₆, cross-check it with G
        @test small_group_identification(G) == (order(G), 1)

        # [[84, 10, 9]]
        l = 7
        m = 3
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^m, r^14, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r^7, r^8, s*r^10]
        b_elts = [one(G), s, r^5, s^2*r^13]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == 2*l*m
        @test describe(G) == "C$l x S$m"
        @test code_n(c) == 84 && code_k(c) == 10
        # Oscar.describe(Oscar.small_group(2*l*m, 3)) is C₇ x S₃, cross-check it with G
        @test small_group_identification(G) == (order(G), 3)

        # [[96, 6, 12]]
        l = 12
        m = 4
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^m, r^l, s^(-1)*r*s*r])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), s, r^9, s * r]
        b_elts = [one(G), s^2 * s^9, r^7, r^2]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l*m
        @test describe(G) == "C$l : C$m"
        @test code_n(c) == 96 && code_k(c) == 6
        # Oscar.describe(Oscar.small_group(l*m, 13)) is C₁₂ x C₄, cross-check it with G
        @test small_group_identification(G) == (order(G), 13)

        # [[96, 12, 10]]
        l = 2
        m = 3
        n = 8
        F = free_group(["r", "s"])
        r, s = gens(F)
        G, = quo(F, [s^(l*m), r^n, r^(-1)*s*r*s])
        F2G = group_algebra(GF(2), G)
        r, s = gens(G)
        a_elts = [one(G), r, s^3 * r^2, s^2 * r^3]
        b_elts = [one(G), r, s^4 * r^6, s^5 * r^3]
        c = get_code(a_elts, b_elts, F2G)
        @test order(G) == l*m*n
        @test describe(G) == "C$l x (C$m : C$n)"
        @test code_n(c) == 96 && code_k(c) == 12
        # Oscar.describe(Oscar.small_group(l*m*n, 9)) is C₂ x (C₃ : C₈), cross-check it with G
        @test small_group_identification(G) == (order(G), 9)
    end
end
