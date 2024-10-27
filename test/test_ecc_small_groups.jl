@testitem "ECC 2BGA Hecke Small Groups" begin
    using Hecke: group_algebra, GF, abelian_group, gens, quo, one, GroupAlgebra, small_group
    using QuantumClifford.ECC
    using QuantumClifford.ECC: code_k, code_n, two_block_group_algebra_codes

    @testset "Hecke Small Groups without extra relations for single cyclic groups" begin
        # [[72, 8, 9]]
        l = 36
        group_id = 2
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^36  ==  1 # presentation ⟨r|r³⁶⟩ satisfied
        a_elts = [one(G), r^28]
        b_elts = [one(G), r, r^18, r^12, r^29, r^14]
        a = sum(GA(x) for x in a_elts)
        b = sum(GA(x) for x in b_elts)
        c = two_block_group_algebra_codes(a,b)
        @test code_n(c) == 72 && code_k(c) == 8

        # [[54, 6, 9]]
        l = 27
        group_id = 1
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^27  ==  1 # presentation ⟨r|r²⁷⟩ satisfied
        a_elts = [one(G), r, r^3, r^7]
        b_elts = [one(G), r, r^12, r^19]
        a = sum(GA(x) for x in a_elts)
        b = sum(GA(x) for x in b_elts)
        c = two_block_group_algebra_codes(a,b)
        @test code_n(c) == 54 && code_k(c) == 6

        # [[60, 6, 10]]
        l = 30
        group_id = 4
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^30  ==  1 # presentation ⟨r|r³⁰⟩ satisfied
        a_elts = [one(G), r^10, r^6, r^13]
        b_elts = [one(G), r^25, r^16, r^12]
        a = sum(GA(x) for x in a_elts)
        b = sum(GA(x) for x in b_elts)
        c = two_block_group_algebra_codes(a,b)
        @test code_n(c) == 60 && code_k(c) == 6

        # [[70, 8, 10]]
        l = 35
        group_id = 1
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^35  ==  1 # presentation ⟨r|r³⁵⟩ satisfied
        a_elts = [one(G), r^15, r^16, r^18]
        b_elts = [one(G), r, r^24, r^27]
        a = sum(GA(x) for x in a_elts)
        b = sum(GA(x) for x in b_elts)
        c = two_block_group_algebra_codes(a,b)
        @test code_n(c) == 70 && code_k(c) == 8

        # [[72, 8, 10]]
        l = 36
        group_id = 2
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^36  ==  1 # presentation ⟨r|r³⁶⟩ satisfied
        a_elts = [one(G), r^9, r^28, r^31]
        b_elts = [one(G), r, r^21, r^34]
        a = sum(GA(x) for x in a_elts)
        b = sum(GA(x) for x in b_elts)
        c = two_block_group_algebra_codes(a,b)
        @test code_n(c) == 72 && code_k(c) == 8

        # [[72, 10, 9]]
        l = 36
        group_id = 2
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^36  ==  1 # presentation ⟨r|r³⁶⟩ satisfied
        a_elts = [one(G), r^9, r^28, r^13]
        b_elts = [one(G), r, r^3, r^22]
        a = sum(GA(x) for x in a_elts)
        b = sum(GA(x) for x in b_elts)
        c = two_block_group_algebra_codes(a,b)
        @test code_n(c) == 72 && code_k(c) == 10
    end
end
