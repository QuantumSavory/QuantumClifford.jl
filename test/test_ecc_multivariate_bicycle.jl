@testitem "ECC Multivaraite Bicycle" begin
    using Hecke
    using JuMP
    using GLPK
    using Hecke: group_algebra, GF, abelian_group, gens
    using QuantumClifford
    using QuantumClifford.ECC: two_block_group_algebra_codes, code_k, code_n, parity_checks

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

    # Multivariate Bicycle codes taken from Table 1 of [voss2024multivariatebicyclecodes](@cite)
    @testset "Weight-4 QLDPC Codes" begin
        # [[112, 8, 5]]
        l=7; m=8
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = z^2 + z^6
        B = x + x^6
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 112 && code_k(c) == 8
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 5

        # [[64, 2, 8]]
        l=8; m=4
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x + x^2
        B = x^3 + y
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 64 && code_k(c) == 2
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 8

        # [[72, 2, 8]]
        l=4; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x + y^2
        B = x^2 + y^2
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 72 && code_k(c) == 2
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 8

        # [[96, 2, 8]]
        l=6; m=8
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^5 + y^6
        B = z + z^4
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 96 && code_k(c) == 2

        # [[112, 2, 10]]
        l=7; m=8
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = z^6 + x^5
        B = z^2 + y^5
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 112 && code_k(c) == 2

        # [[144, 2, 12]]
        l=8; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^3 + y^7
        B = x + y^5
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 144 && code_k(c) == 2
    end

    @testset "Weight-5 QLDPC Codes" begin
        # [[30, 4, 5]]
        l=3; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x + z^4
        B = x + y^2 + z^2 
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 30 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 5

        # [[72, 4, 8]]
        l=4; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x + y^3
        B = x^2 + y + y^2 
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 72 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 8

        # [96, 4, 8]]
        l=8; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^6 + x^3
        B = z^5 + x^5 + y 
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 96 && code_k(c) == 4
    end

    @testset "Weight-6 QLDPC codes" begin
        # [[30, 6, 4]]
        l=5; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^4 + z^3
        B = x^4 + x + z^4 + y
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 30 && code_k(c) == 6
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 4

        # [[48, 6, 6]]
        l=4; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^2 + y^4
        B = x^3 + z^3 + y^2 + y
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 48 && code_k(c) == 6
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 6

        # [[40, 4, 6]]
        l=4; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^2 + y
        B = y^4 + y^2 + x^3 + x
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 40 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 6

        # [[48, 4, 6]]
        l=4; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^3 + y^5
        B = x + z^5 + y^5 + y^2
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 48 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 6
    end

    @testset "Weight-7 QLDPC codes" begin
        # [[30, 4, 5]]
        l=5; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^4 + x^2
        B = x + x^2 + y + z^2 + z^3
        c = two_block_group_algebra_codes(A, B)
        @test code_n(c) == 30 && code_k(c) == 4
        hx, lx = get_hx_lx(c)
        @test code_distance(hx, lx, code_k(c)) == 5
    end
end
