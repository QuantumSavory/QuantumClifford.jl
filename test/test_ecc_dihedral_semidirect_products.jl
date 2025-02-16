@testitem "ECC Dihedral group via semidirect_product" begin
    using Nemo: FqFieldElem
    using Hecke: group_algebra, GF, abelian_group, gens, quo, one, GroupAlgebra, GroupAlgebraElem
    using QuantumClifford.ECC
    using QuantumClifford.ECC: code_k, code_n, two_block_group_algebra_codes
    using Oscar: small_group_identification, describe, order, FPGroupElem, FPGroup, FPGroupElem, semidirect_product, automorphism_group, hom, gen, cyclic_group, SemidirectProductGroup, PcGroup, BasicGAPGroupElem, normal_subgroup

    function get_code(a_elts::Vector{GroupAlgebraElem{FqFieldElem, GroupAlgebra{FqFieldElem, SemidirectProductGroup{PcGroup, PcGroup}, BasicGAPGroupElem{SemidirectProductGroup{PcGroup, PcGroup}}}}}, b_elts::Vector{GroupAlgebraElem{FqFieldElem, GroupAlgebra{FqFieldElem, SemidirectProductGroup{PcGroup, PcGroup}, BasicGAPGroupElem{SemidirectProductGroup{PcGroup, PcGroup}}}}}, GA::GroupAlgebra{FqFieldElem, SemidirectProductGroup{PcGroup, PcGroup}, BasicGAPGroupElem{SemidirectProductGroup{PcGroup, PcGroup}}})
        a = sum(GA(x) for x in a_elts)
        b = sum(GA(x) for x in b_elts)
        c = two_block_group_algebra_codes(a,b)
        return c
    end

    function semidirectproduct(m::Int)
        Cₘ	= cyclic_group(m)
        C₂	= cyclic_group(2)
        A	= automorphism_group(Cₘ)
        # Given specific Dihedral group presentation, choose r -> r⁻¹
        au	= A(hom(Cₘ,Cₘ,[Cₘ[1]],[Cₘ[1]^-1]))
        f	= hom(C₂,A,[C₂[1]],[au])
        G 	= semidirect_product(Cₘ,f,C₂)
        s	= gen(G, 1)
        r	= gen(G, 2)
        @test	r^m == s^2 == (r*s)^2
        GA	= group_algebra(GF(2), G)
        r, s	= gens(GA)[2], gens(GA)[3];
        return r, s, GA, G
    end

    @testset "Reproduce Table 3 of lin2024quantum" begin
        # [[24, 8, 3]]
        m = 6
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^4]
        b_elts = [one(r), s*r^4, r^3, r^4, s*r^2, r]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test describe(normal_subgroup(G)) == "C$m"
        @test code_n(c) == 24 && code_k(c) == 8
        # Oscar.describe(Oscar.small_group(2*m, 4)) is D₁₂, cross-check it with G
        @test small_group_identification(G) == (order(G), 4)

        # [[24, 12, 2]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^3]
        b_elts = [one(r), s*r, r^3, r^4, s*r^4, r]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 24 && code_k(c) == 12
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 4)) is D₁₂, cross-check it with G
        @test small_group_identification(G) == (order(G), 4)

        # [[32, 8, 4]]
        m = 8
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^2]
        b_elts = [one(r), s*r^5, s*r^4, r^2, s*r^7, s*r^6]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 32 && code_k(c) == 8
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 7)) is D₁₆, cross-check it with G
        @test small_group_identification(G) == (order(G), 7)

        # [[32, 16, 2]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^4]
        b_elts = [one(r), s*r^3, s*r^6, r^4, s*r^7, s*r^2]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 32 && code_k(c) == 16
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 7)) is D₁₆, cross-check it with G
        @test small_group_identification(G) == (order(G), 7)

        # [[36, 12, 3]]
        m = 9
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^3]
        b_elts = [one(r), s, r, r^3, s*r^3, r^4]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 36 && code_k(c) == 12
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 1)) is D₁₈, cross-check it with G
        @test small_group_identification(G) == (order(G), 1)

        # [[40, 8, 5]]
        m = 10
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^2]
        b_elts = [one(r), s*r^4, r^5, r^2, s*r^6, r]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 40 && code_k(c) == 8
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 4)) is D₂₀, cross-check it with G
        @test small_group_identification(G) == (order(G), 4)

        # [[40, 20, 2]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^5]
        b_elts = [one(r), s*r^2, r^5, r^6, s*r^7, r]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 40 && code_k(c) == 20
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 4)) is D₂₀, cross-check it with G
        @test small_group_identification(G) == (order(G), 4)

        # [[48, 8, 6]]
        m =  12
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^10]
        b_elts = [one(r), s*r^8, r^9, r^4, s*r^2, r^5]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 48 && code_k(c) == 8
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 6)) is D₂₄, cross-check it with G
        @test small_group_identification(G) == (order(G), 6)

        # [[48, 12, 4]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^3]
        b_elts = [one(r), s*r^7, r^3, r^4, s*r^10, r^7]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 48 && code_k(c) == 12
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 6)) is D₂₄, cross-check it with G
        @test small_group_identification(G) == (order(G), 6)

        # [[48, 16, 3]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^8]
        b_elts = [one(r), s*r^8, r^9, r^8, s*r^4, r^5]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 48 && code_k(c) == 16
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 6)) is D₂₄, cross-check it with G
        @test small_group_identification(G) == (order(G), 6)

        # [[48, 24, 2]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^6]
        b_elts = [one(r), s*r^11, r^6, s*r^5, r, r^7]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 48 && code_k(c) == 24
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 6)) is D₂₄, cross-check it with G
        @test small_group_identification(G) == (order(G), 6)

        # [[56, 8, 7]]
        m = 14
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^4]
        b_elts = [one(r), s*r^11, r^7, s*r^5, r^12, r^9]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 56 && code_k(c) == 8
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 3)) is D₂₈, cross-check it with G
        @test small_group_identification(G) == (order(G), 3)

        # [[56, 28, 2]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^7]
        b_elts = [one(r), s*r^2, r^7, r^8, s*r^9, r]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 56 && code_k(c) == 28
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 3)) is D₂₈, cross-check it with G
        @test small_group_identification(G) == (order(G), 3)

        # [[60, 12, 5]]
        m = 15
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^12]
        b_elts = [one(r), s*r^14, r^5, r^12, s*r^11, r^14]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 60 && code_k(c) == 12
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 3)) is D₃₀, cross-check it with G
        @test small_group_identification(G) == (order(G), 3)

        # [[60, 20, 3]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^5]
        b_elts = [one(r), s*r^13, r^5, r^12, s*r^3, r^2]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 60 && code_k(c) == 20
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 3)) is D₃₀, cross-check it with G
        @test small_group_identification(G) == (order(G), 3)

        # [[64, 8, 8]]
        m = 16
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^6]
        b_elts = [one(r), s*r^12, s*r^9, r^6, s, s*r]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 64 && code_k(c) == 8
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 18)) is D₃₂, cross-check it with G
        @test small_group_identification(G) == (order(G), 18)

        # [[64, 16, 8]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^4]
        b_elts = [one(r), s*r^10, s*r^3, r^4, s*r^14, s*r^7]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 64 && code_k(c) == 16
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 18)) is D₃₂, cross-check it with G
        @test small_group_identification(G) == (order(G), 18)

        # [[64, 32, 2]]
        r, s, GA, G = semidirectproduct(m)
        a_elts = [one(r), r^8]
        b_elts = [one(r), s*r^11, s*r^12, r^8, s*r^3, s*r^4]
        c = get_code(a_elts, b_elts, GA)
        @test order(G) == 2*m
        @test describe(G) == "D$(m*2)"
        @test code_n(c) == 64 && code_k(c) == 32
        @test describe(normal_subgroup(G)) == "C$m"
        # Oscar.describe(Oscar.small_group(2*m, 18)) is D₃₂, cross-check it with G
        @test small_group_identification(G) == (order(G), 18)
    end
end
