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
        A = 1 + r^28
        B = 1 + r + r^18 + r^12 + r^29 + r^14
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 72 && code_k(c) == 8

        # [[54, 6, 9]]
        l = 27
        group_id = 1
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^27  ==  1 # presentation ⟨r|r²⁷⟩ satisfied
        A = 1 + r + r^3  + r^7
        B = 1 + r + r^12 + r^19
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 54 && code_k(c) == 6

        # [[60, 6, 10]]
        l = 30
        group_id = 4
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^30  ==  1 # presentation ⟨r|r³⁰⟩ satisfied
        A = 1 + r^10 + r^6  + r^13
        B = 1 + r^25 + r^16 + r^12
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 60 && code_k(c) == 6

        # [[70, 8, 10]]
        l = 35
        group_id = 1
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^35  ==  1 # presentation ⟨r|r³⁵⟩ satisfied
        A = 1 + r^15 + r^16 + r^18
        B = 1 + r    + r^24 + r^27
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 70 && code_k(c) == 8

        # [[72, 8, 10]]
        l = 36
        group_id = 2
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^36  ==  1 # presentation ⟨r|r³⁶⟩ satisfied
        A = 1 + r^9 + r^28 + r^31
        B = 1 + r   + r^21 + r^34
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 72 && code_k(c) == 8

        # [[72, 10, 9]]
        l = 36
        group_id = 2
        G = small_group(l, group_id)
        GA = group_algebra(GF(2), G)
        r = prod(gens(GA))
        @test r^36  ==  1 # presentation ⟨r|r³⁶⟩ satisfied
        A = 1 + r^9 + r^28 + r^13
        B = 1 + r   + r^3  + r^22
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 72 && code_k(c) == 10
    end
end
