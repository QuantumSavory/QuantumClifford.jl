@testitem "ECC coprime Bivaraite Bicycle" tags=[:ecc, :ecc_bespoke_checks] begin
    using Nemo
    using Nemo: gcd
    using Hecke
    using JuMP
    using HiGHS
    using Hecke: group_algebra, GF, abelian_group, gens
    using QuantumClifford.ECC: two_block_group_algebra_code, code_k, code_n, distance, DistanceMIPAlgorithm

    @testset "Reproduce Table 2 wang2024coprime" begin
        # [[30,4,6]]
        l=3; m=5;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋   + 𝜋^2
        B = 𝜋 + 𝜋^3 + 𝜋^8
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6
        @test code_n(c) == 30 && code_k(c) == 4

        # [[42,6,6]]
        l=3; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋^2 + 𝜋^3
        B = 𝜋 + 𝜋^3 + 𝜋^11
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6
        @test code_n(c) == 42 && code_k(c) == 6

        # [[70,6,8]]
        l=5; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^5;
        B = 1 + 𝜋 + 𝜋^12;
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8
        @test code_n(c) == 70 && code_k(c) == 6

        # [[108,12,6]]
        l=2; m=27;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 𝜋^2 + 𝜋^5  + 𝜋^44
        B = 𝜋^8 + 𝜋^14 + 𝜋^47
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        i = rand(1:code_k(c))
        @test distance(c, DistanceMIPAlgorithm(logical_qubit=i; solver=HiGHS)) == 6
        @test code_n(c) == 108 && code_k(c) == 12

        # [[126,12,10]]
        l=7; m=9
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1   + 𝜋    + 𝜋^58
        B = 𝜋^3 + 𝜋^16 + 𝜋^44
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 126 && code_k(c) == 12

        # [[154,6,16]]
        l=7; m=11
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^31
        B = 1 + 𝜋^19 + 𝜋^53
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 154 && code_k(c) == 6
    end

    @testset "Reproduce Appendix B Table 4 wang2024coprime" begin
        # Weight-6 coprime-BB codes from Appendix B, Table 4

        # [[28,6,4]]
        l=2; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3
        B = 1 + 𝜋 + 𝜋^10
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 28 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

        # [[36,8,4]]
        l=2; m=9;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋^2 + 𝜋^10
        B = 1 + 𝜋^4 + 𝜋^8
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 36 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

        # [[42,10,4]]
        l=3; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^5
        B = 1 + 𝜋^2 + 𝜋^10
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 42 && code_k(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

        # [[48,4,8]]
        l=3; m=8;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^2
        B = 1 + 𝜋^2 + 𝜋^10
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 48 && code_k(c) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8

        # [[56,6,8]]
        l=4; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3
        B = 1 + 𝜋^5 + 𝜋^11
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 56 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8

        # [[60,16,4]]
        l=3; m=10;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋^2 + 𝜋^8
        B = 1 + 𝜋^4 + 𝜋^16
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 60 && code_k(c) == 16
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

        # [[66,4,10]]
        l=3; m=11;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^5
        B = 1 + 𝜋 + 𝜋^23
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 66 && code_k(c) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 10

        # [[84,6,10]]
        l=6; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3
        B = 1 + 𝜋^8 + 𝜋^31
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 84 && code_k(c) == 6
        i = rand(1:code_k(c))
        @test distance(c, DistanceMIPAlgorithm(logical_qubit=i; solver=HiGHS)) == 10

        # NOTE: [[90,4,12]] from Table 4 with l=5, m=9, a(π)=1+π+π⁴, b(π)=1+π⁸+π³⁴
        # is not included because it produces k=0 instead of the paper's claimed k=4.

        # [[90,8,8]]
        l=5; m=9;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^12
        B = 1 + 𝜋^2 + 𝜋^9
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 90 && code_k(c) == 8
        i = rand(1:code_k(c))
        @test distance(c, DistanceMIPAlgorithm(logical_qubit=i; solver=HiGHS)) == 8

        # [[112,6,12]]
        l=7; m=8;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3
        B = 1 + 𝜋^5 + 𝜋^25
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 112 && code_k(c) == 6

        # [[126,6,14]]
        l=7; m=9;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋^4 + 𝜋^19
        B = 1 + 𝜋^6 + 𝜋^16
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 126 && code_k(c) == 6

        # [[132,4,14]]
        l=6; m=11;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^2
        B = 1 + 𝜋^11 + 𝜋^28
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 132 && code_k(c) == 4

        # [[180,8,16]]
        l=9; m=10;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^4
        B = 1 + 𝜋^23 + 𝜋^62
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 180 && code_k(c) == 8
    end

    @testset "Reproduce Appendix B Table 5 wang2024coprime" begin
        # Weight-8 coprime-BB codes from Appendix B, Table 5

        # [[24,8,4]]
        l=3; m=4;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3 + 𝜋^4
        B = 1 + 𝜋^2 + 𝜋^5 + 𝜋^9
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 24 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

        # [[30,10,4]]
        l=3; m=5;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^2 + 𝜋^7
        B = 1 + 𝜋 + 𝜋^4 + 𝜋^10
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 30 && code_k(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

        # [[30,6,5]]
        l=3; m=5;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3 + 𝜋^4
        B = 1 + 𝜋 + 𝜋^3 + 𝜋^7
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 30 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5

        # [[40,14,4]]
        l=4; m=5;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^6 + 𝜋^15
        B = 1 + 𝜋^2 + 𝜋^5 + 𝜋^7
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 40 && code_k(c) == 14
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

        # [[40,6,6]]
        l=4; m=5;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^2 + 𝜋^3
        B = 1 + 𝜋 + 𝜋^3 + 𝜋^10
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 40 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

        # [[40,8,5]]
        l=4; m=5;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^4 + 𝜋^5
        B = 1 + 𝜋 + 𝜋^4 + 𝜋^9
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 40 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5

        # [[42,12,5]]
        l=3; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3 + 𝜋^13
        B = 1 + 𝜋 + 𝜋^4 + 𝜋^9
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 42 && code_k(c) == 12
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5

        # [[42,6,7]]
        l=3; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3 + 𝜋^4
        B = 1 + 𝜋 + 𝜋^6 + 𝜋^10
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 42 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 7

        # [[48,6,8]]
        l=3; m=8;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^2 + 𝜋^3
        B = 1 + 𝜋^3 + 𝜋^9 + 𝜋^14
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 48 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8

        # [[48,10,6]]
        l=3; m=8;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3 + 𝜋^10
        B = 1 + 𝜋^3 + 𝜋^7 + 𝜋^16
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 48 && code_k(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

        # [[56,8,8]]
        l=4; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^2 + 𝜋^4
        B = 1 + 𝜋^2 + 𝜋^6 + 𝜋^19
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 56 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8

        # [[56,14,6]]
        l=4; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^4 + 𝜋^9
        B = 1 + 𝜋 + 𝜋^17 + 𝜋^20
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 56 && code_k(c) == 14
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

        # [[60,12,7]]
        l=5; m=6;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^2 + 𝜋^7
        B = 1 + 𝜋^3 + 𝜋^12 + 𝜋^25
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 60 && code_k(c) == 12
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 7

        # [[60,6,9]]
        l=5; m=6;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^3 + 𝜋^4
        B = 1 + 𝜋^2 + 𝜋^11 + 𝜋^18
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 60 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 9

        # [[70,8,9]]
        l=5; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^2 + 𝜋^4
        B = 1 + 𝜋 + 𝜋^6 + 𝜋^24
        c = two_block_group_algebra_code(A, B)
        @test gcd([l,m]) == 1
        @test code_n(c) == 70 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 9
    end
end
