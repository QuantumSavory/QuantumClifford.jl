@testitem "ECC QLDPC minimum distance" begin
    using Hecke
    using JuMP
    using HiGHS
    using Hecke: group_algebra, GF, abelian_group, gens
    using QuantumClifford.ECC: two_block_group_algebra_codes, generalized_bicycle_codes, code_k, code_n, distance, DistanceMIPAlgorithm

    @testset "minimum distance properties: 2BGA" begin
        # [[56, 8, 7]] 2BGA code code taken from Table 2 of [lin2024quantum](@cite)
        # m = 14
        GA = group_algebra(GF(2), abelian_group([14,2]))
        x, s = gens(GA)
        A = 1 + x^8
        B = 1 + x^7 + s + x^8 + x^9 + s*x^4
        c = two_block_group_algebra_codes(A,B)
        # [[56, 8, 7]] 2BGA code
        # minimum distance is exact, d = 7
        i = rand(1:code_k(c))
        # By default, the minimum distance for the Z-type logical operator is computed.
        # The minimum distance for X-type logical operators is the same.
        @test distance(c, DistanceMIPAlgorithm(logical_qubit=i; solver=HiGHS)) == distance(c, DistanceMIPAlgorithm(logical_qubit=i, logical_operator_type=:Z; solver=HiGHS)) == 7
    end

    @testset "minimum distance properties: GB" begin
        # [48, 6, 8]] GB code, # minimum distance is exact, d = 8
        l = 24
        c = generalized_bicycle_codes([0, 2, 8, 15], [0, 2, 12, 17], l)
        # minimum distance is exact, d = 8
        i = rand(1:code_k(c))
        # By default, the minimum distance for the Z-type logical operator is computed.
        # The minimum distance for X-type logical operators is the same.
        @test distance(c, DistanceMIPAlgorithm(logical_qubit=i; solver=HiGHS)) == distance(c, DistanceMIPAlgorithm(logical_qubit=i, logical_operator_type=:Z; solver=HiGHS)) == 8
    end

    @testset "minimum distance properties: BB" begin
        # [[72, 12, 6]] BB code from Table 3 [bravyi2024high](@cite)
        l=6; m=6
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        A = x^3 + y + y^2
        B = y^3 + x + x^2
        c = two_block_group_algebra_codes(A,B)
        # minimum distance is exact, d = 6
        i = rand(1:code_k(c))
        # By default, the minimum distance for the Z-type logical operator is computed.
        # The minimum distance for X-type logical operators is the same.
        @test distance(c, DistanceMIPAlgorithm(logical_qubit=i; solver=HiGHS)) == distance(c, DistanceMIPAlgorithm(logical_qubit=i, logical_operator_type=:Z; solver=HiGHS)) == 6
    end

    @testset "minimum distance properties: coprime BB" begin
        # [[70, 6, 8]] coprime BB code from Table 2 [wang2024coprime](@cite)
        l=5; m=7;
        GA = group_algebra(GF(2), abelian_group([l*m]))
        𝜋 = gens(GA)[1]
        A = 1 + 𝜋 + 𝜋^5;
        B = 1 + 𝜋 + 𝜋^12;
        c = two_block_group_algebra_codes(A, B)
        # minimum distance is exact, d = 8
        i = rand(1:code_k(c))
        @test code_n(c) == 70 && code_k(c) == 6
        # By default, the minimum distance for the Z-type logical operator is computed.
        # The minimum distance for X-type logical operators is the same.
        @test distance(c, DistanceMIPAlgorithm(logical_qubit=i; solver=HiGHS)) == distance(c, DistanceMIPAlgorithm(logical_qubit=i, logical_operator_type=:Z; solver=HiGHS)) == 8
    end

    @testset "minimum distance properties: Weight-7 MB" begin
        # [[30, 4, 5]] MB code from Table 1 of [voss2024multivariatebicyclecodes](@cite)
        l=5; m=3
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        z = x*y
        A = x^4 + x^2
        B = x + x^2 + y + z^2 + z^3
        c = two_block_group_algebra_codes(A, B)
        # minimum distance is exact, d = 5
        i = rand(1:code_k(c))
        @test code_n(c) == 30 && code_k(c) == 4
        # By default, the minimum distance for the Z-type logical operator is computed.
        # The minimum distance for X-type logical operators is the same.
        @test distance(c, DistanceMIPAlgorithm(logical_qubit=i; solver=HiGHS)) == distance(c, DistanceMIPAlgorithm(logical_qubit=i, logical_operator_type=:Z; solver=HiGHS)) == 5
    end
end
