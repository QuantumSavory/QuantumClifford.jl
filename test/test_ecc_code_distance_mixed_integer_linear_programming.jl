@testitem "ECC Minimum Distance Solver MILP" begin
    using Hecke
    using JuMP
    using GLPK
    using Hecke: group_algebra, GF, abelian_group, gens, one
    using QuantumClifford
    using QuantumClifford.ECC: two_block_group_algebra_codes, code_k, code_n, parity_checks

    include("test_ecc_util.jl") # minimum_distance

    function get_hx_lx(c)
        hx = stab_to_gf2(Stabilizer(parity_checks(c))[1:end÷2,:])
        lx = stab_to_gf2(logicalxview(canonicalize!(MixedDestabilizer(parity_checks(c)))))
        return hx, lx
    end

    @testset "Reproduce Table 3 bravyi2024high" begin
        # [[72, 12, 6]]
        l=6; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 6

        # [[90, 8, 10]]
        l=15; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^9 + y   + y^2
        B = 1   + x^2 + x^7
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 10

        # [[108, 8, 10]]
        l=9; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 10
    end

    @testset "Reproduce Table 1 berthusen2024toward" begin
        # [[72, 8, 6]]
        l=12; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^9 + y + y^2
        B = 1   + x + x^11
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 6

        # [[90, 8, 6]]
        l=9; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^8 + y^4 + y
        B = y^5 + x^8 + x^7
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 6

        # [[120, 8, 8]]
        l=12; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^10 + y^4 + y
        B = 1    + x   + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 8
    end

    @testset "Reproduce Table 1 wang2024coprime" begin
        # [[54, 8, 6]]
        l=3; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y^2 + y^4
        B = y^3 + x   + x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 6

        # [[98, 6, 12]]
        l=7; m=7
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y^5 + y^6
        B = y^2 + x^3 + x^5
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 12

        # [[126, 8, 10]]
        l=3; m=21
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y^2 + y^10
        B = y^3 + x  +  x^2
        c = two_block_group_algebra_codes(A,B)
        hx, lx = get_hx_lx(c)
        i = rand(1:code_k(c))
        @test minimum_distance(hx, lx[i, :]) == 10
    end
end
