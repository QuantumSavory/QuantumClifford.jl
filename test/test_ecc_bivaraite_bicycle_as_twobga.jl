@testitem "ECC Bivaraite Bicycle as 2BGA" tags=[:ecc] begin
    using Hecke
    using HiGHS
    using JuMP
    using Hecke: group_algebra, GF, abelian_group, gens, one
    using QuantumClifford.ECC: two_block_group_algebra_codes, code_k, code_n, distance, DistanceMIPAlgorithm

    @testset "Reproduce Table 3 bravyi2024high" begin
        # [[72, 12, 6]]
        l=6; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 72 && code_k(c) == 12
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

        # [[90, 8, 10]]
        l=15; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^9 + y   + y^2
        B = 1   + x^2 + x^7
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 90 && code_k(c) == 8

        # [[108, 8, 10]]
        l=9; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 108 && code_k(c) == 8

        # [[144, 12, 12]]
        l=12; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 144 && code_k(c) == 12

        # [[288, 12, 12]]
        l=12; m=12
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y^2 + y^7
        B = y^3 + x   + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 288 && code_k(c) == 12

        # [[360, 12, ≤ 24]]
        l=30; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^9 + y    + y^2
        B = y^3 + x^25 + x^26
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 360 && code_k(c) == 12

        # [[756, 16, ≤ 34]]
        l=21; m=18
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y^10 + y^17
        B = y^5 + x^3  + x^19
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 756 && code_k(c) == 16
    end

    @testset "Reproduce Table 1 berthusen2024toward" begin
        # [[72, 8, 6]]
        l=12; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^9 + y + y^2
        B = 1   + x + x^11
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 72 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

        # [[90, 8, 6]]
        l=9; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^8 + y^4 + y
        B = y^5 + x^8 + x^7
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 90 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

        # [[120, 8, 8]]
        l=12; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^10 + y^4 + y
        B = 1    + x   + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 120 && code_k(c) == 8

        # [[150, 8, 8]]
        l=15; m=5
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^5 + y^2 + y^3
        B = y^2 + x^7 + x^6
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 150 && code_k(c) == 8

        # [[196, 12, 8]]
        l=14; m=7
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^6 + y^5 + y^6
        B = 1   + x^4 + x^13
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 196 && code_k(c) == 12
    end

    @testset "Reproduce Table 1 wang2024coprime" begin
        # [[54, 8, 6]]
        l=3; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y^2 + y^4
        B = y^3 + x   + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 54 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

        # [[98, 6, 12]]
        l=7; m=7
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y^5 + y^6
        B = y^2 + x^3 + x^5
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 98 && code_k(c) == 6

        # [[126, 8, 10]]
        l=3; m=21
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y^2 + y^10
        B = y^3 + x  +  x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 126 && code_k(c) == 8

        # [[150, 16, 8]]
        l=5; m=15
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y^6 + y^8
        B = y^5 + x   + x^4
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 150 && code_k(c) == 16

        # [[162, 8, 14]]
        l=3; m=27
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1    + y^10 + y^14
        B = y^12 + x    + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 162 && code_k(c) == 8

        # [[180, 8, 16]]
        l=6; m=15
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y   + y^2
        B = y^6 + x^4 + x^5
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 180 && code_k(c) == 8
    end

    @testset "Reproduce Table 1 eberhardt2024logical" begin
        # [[108, 16, 6]]
        l=6; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1 +   y   + y^2
        B = y^3 + x^2 + x^4
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 108 && code_k(c) == 16

        # [[128, 14, 12]]
        l=8; m=8
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^2 + y + y^3 + y^4
        B = y^2 + x + x^3 + x^4
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 128 && code_k(c) == 14

        # [[162, 4, 16]]
        l=9; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + x + y
        B = x^3 + y + y^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 162 && code_k(c) == 4

        # [[162, 12, 8]]
        l=9; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + x   + y^6
        B = y^3 + x^2 + x^3
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 162 && code_k(c) == 12

        # [[162, 24, 6]]
        l=9; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = 1   + y   + y^2
        B = y^3 + x^3 + x^6
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 162 && code_k(c) == 24

        # [[270, 8, 18]]
        l=9; m=15
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 270 && code_k(c) == 8

        # [[98, 6, 12]]
        l=7; m=7
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x + y^3 + y^4
        B = y + x^3 + x^4
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 98 && code_k(c) == 6

        # [[162, 8, 12]]
        l=9; m=9
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        @test code_n(c) == 162 && code_k(c) == 8
    end
end
