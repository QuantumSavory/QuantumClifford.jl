@testitem "ECC 2BGA Hecke Small Groups" begin
    @static if Sys.iswindows() && VERSION < v"1.11"
        import Hecke: group_algebra, GF, abelian_group, gens, quo, one, small_group, prod
        using QuantumClifford.ECC
        using QuantumClifford.ECC: code_k, code_n, two_block_group_algebra_codes

        @testset "Hecke Small Groups for single cyclic groups" begin
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

        @testset "Hecke Small Groups Block 4 Table I of arxiv:2306.16400" begin
            # [[72, 8, 9]]
            m = 4
            n = 9
            l = 36
            group_id = 1
            G = small_group(l, group_id)
            GA = group_algebra(GF(2), G)
            s, r = gens(GA)[1], gens(GA)[2]
            @test s^m == r^n == s^-1*r*s*r
            a = 1 + s + r + s*r^6
            b = 1 + s^2*r + s^2*r^6 + r^2
            c = two_block_group_algebra_codes(a,b)
            @test code_n(c) == 72 && code_k(c) == 8

            # [[80, 8, 10]]
            m = 5
            n = 8
            l = 40
            group_id = 1
            G = small_group(l, group_id)
            GA = group_algebra(GF(2), G)
            s, r = gens(GA)[2], gens(GA)[1]
            @test s^m == r^n == r^-1*s*r*s
            a = 1 + r + s + s^3*r^5
            b = 1 + r^2 + s*r^4 + s^3*r^2
            c = two_block_group_algebra_codes(a,b)
            @test code_n(c) == 80 && code_k(c) == 8

            # [[96, 10, 12]]
            m = 4
            n = 12
            l = 48
            group_id = 13
            G = small_group(l, group_id)
            GA = group_algebra(GF(2), G)
            s = gens(GA)[1]
            r = gens(GA)[2]*gens(GA)[3]
            @test s^m == r^n == s^-1*r*s*r
            a = 1 + s + r^9 + s*r
            b = 1 + s^2*r^9 + r^7 + r^2
            c = two_block_group_algebra_codes(a,b)
            @test code_n(c) == 96 && code_k(c) == 10

            # [[96, 11, 9]]
            m = 2
            n = 24
            l = 48
            group_id = 5
            G = small_group(l, group_id)
            GA = group_algebra(GF(2), G)
            r = gens(GA)[2]*gens(GA)[3]
            s = gens(GA)[1];
            @test s^m == r^n == (r*s)^8
            a = 1 + s + r^9 + s*r^13
            b = 1 + r^9 + s*r^18 + r^7
            c = two_block_group_algebra_codes(a,b)
            @test code_n(c) == 96 && code_k(c) == 11

            # [[96, 12, 10]]
            m = 6
            n = 8
            l = 48
            group_id = 9
            G = small_group(l, group_id)
            GA = group_algebra(GF(2), G)
            r = gens(GA)[1]*gens(GA)[2]
            s = gens(GA)[3]
            @test s^m == r^n == r^-1*s*r*s
            a = 1 + r + s^3*r^2 + s^2*r^3
            b = 1 + r + s^4*r^6 + s^5*r^3
            c = two_block_group_algebra_codes(a,b)
            @test code_n(c) == 96 && code_k(c) == 12
        end
    end
end
