@testitem "ECC 2BGA" begin
    using Hecke
    using JuMP
    using GLPK
    using Hecke: group_algebra, GF, abelian_group, gens, one
    using QuantumClifford
    using QuantumClifford.ECC: two_block_group_algebra_codes, code_k, code_n, parity_checks, parity_checks_x

    include("test_ecc_util.jl") # minimum_distance

    function get_hx_lx(c)
        hx = stab_to_gf2(Stabilizer(parity_checks(c))[1:end÷2,:])
        lx = stab_to_gf2(logicalxview(canonicalize!(MixedDestabilizer(parity_checks(c)))))
        return hx, lx
    end

    function code_distance(hx, lx, k)
        w_values = []
        for i in 1:k
            w = minimum_distance(hx, lx[i, :])
            push!(w_values, w)
        end
        return round(Int, sum(w_values)/k)
    end

    @testset "Reproduce Table 2 lin2024quantum" begin
        # codes taken from Table 2 of [lin2024quantum](@cite)

        # m = 4
        GA = group_algebra(GF(2), abelian_group([4,2]))
        x, s = gens(GA)
        A = 1 + x
        B = 1 + x + s + x^2 + s*x + s*x^3
        c = two_block_group_algebra_codes(A,B)
        # [[16, 2, 4]] 2BGA code
        @test code_n(c) == 16 && code_k(c) == 2
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 4
        A = 1 + x
        B = 1 + x + s + x^2 + s*x + x^3
        c = two_block_group_algebra_codes(A,B)
        # [[16, 4, 4]] 2BGA code
        @test code_n(c) == 16 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 4
        A = 1 + s
        B = 1 + x + s + x^2 + s*x + x^2
        c = two_block_group_algebra_codes(A,B)
        # [[16, 8, 2]] 2BGA code
        @test code_n(c) == 16 && code_k(c) == 8
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 2

        # m = 6
        GA = group_algebra(GF(2), abelian_group([6,2]))
        x, s = gens(GA)
        A = 1 + x
        B = 1 + x^3 + s + x^4 + x^2 + s*x
        c = two_block_group_algebra_codes(A,B)
        # [[24, 4, 5]] 2BGA code
        @test code_n(c) == 24 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 5
        A = 1 + x^3
        B = 1 + x^3 + s + x^4 + s*x^3 + x
        c = two_block_group_algebra_codes(A,B)
        # [[24, 12, 2]] 2BGA code
        @test code_n(c) == 24 && code_k(c) == 12
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 2

        # m = 8
        GA = group_algebra(GF(2), abelian_group([8,2]))
        x, s = gens(GA)
        A = 1 + x^6
        B = 1 + s*x^7 + s*x^4 + x^6 + s*x^5 + s*x^2
        c = two_block_group_algebra_codes(A,B)
        # [[32, 8, 4]] 2BGA code
        @test code_n(c) == 32 && code_k(c) == 8
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 4
        A = 1 + s*x^4
        B = 1 + s*x^7 + s*x^4 + x^6 + x^3 + s*x^2
        c = two_block_group_algebra_codes(A,B)
        # [[32, 16, 2]] 2BGA code
        @test code_n(c) == 32 && code_k(c) == 16
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 2

        # m = 10
        GA = group_algebra(GF(2), abelian_group([10,2]))
        x, s = gens(GA)
        A = 1 + x
        B = 1 + x^5 + x^6 + s*x^6 + x^7 + s*x^3
        c = two_block_group_algebra_codes(A,B)
        # [[40, 4, 8]] 2BGA code
        @test code_n(c) == 40 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 8
        A = 1 + x^6
        B = 1 + x^5 + s + x^6 + x + s*x^2
        c = two_block_group_algebra_codes(A,B)
        # [[40, 8, 5]] 2BGA code
        @test code_n(c) == 40 && code_k(c) == 8
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 5
        A = 1 + x^5
        B = 1 + x^5 + s + x^6 + s*x^5 + x
        c = two_block_group_algebra_codes(A,B)
        # [[40, 20, 2]] 2BGA code
        @test code_n(c) == 40 && code_k(c) == 20
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 2

        # m = 12
        GA = group_algebra(GF(2), abelian_group([12,2]))
        x, s = gens(GA)
        A = 1 + s*x^10
        B = 1 + x^3 + s*x^6 + x^4 + x^7 + x^8
        c = two_block_group_algebra_codes(A,B)
        # [[48, 8, 6]] 2BGA code
        @test code_n(c) == 48 && code_k(c) == 8
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 6
        A = 1 + x^3
        B = 1 + x^3 + s*x^6 + x^4 + s*x^9 + x^7
        c = two_block_group_algebra_codes(A,B)
        # [[48, 12, 4]] 2BGA code
        @test code_n(c) == 48 && code_k(c) == 12
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 4
        A = 1 + x^4
        B = 1 + x^3 + s*x^6 + x^4 + x^7 + s*x^10
        c = two_block_group_algebra_codes(A,B)
        # [[48, 16, 3]] 2BGA code
        @test code_n(c) == 48 && code_k(c) == 16
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 3
        A = 1 + s*x^6
        B = 1 + x^3 + s*x^6 + x^4 + s*x^9 + s*x^10
        c = two_block_group_algebra_codes(A,B)
        # [[48, 24, 2]] 2BGA code
        @test code_n(c) == 48 && code_k(c) == 24
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 2

        # m = 14
        GA = group_algebra(GF(2), abelian_group([14,2]))
        x, s = gens(GA)
        A = 1 + x^8
        B = 1 + x^7 + s + x^8 + x^9 + s*x^4
        c = two_block_group_algebra_codes(A,B)
        # [[56, 8, 7]] 2BGA code
        @test code_n(c) == 56 && code_k(c) == 8
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 7
        A = 1 + x^7
        B = 1 + x^7 + s + x^8 + s*x^7 + x
        c = two_block_group_algebra_codes(A,B)
        # [[56, 28, 2]] 2BGA code
        @test code_n(c) == 56 && code_k(c) == 28
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 2
    end
end
