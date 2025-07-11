@testitem "ECC QLDPC minimum distance" tags=[:ecc] begin
    using Hecke
    using JuMP
    using HiGHS
    using Hecke: group_algebra, GF, abelian_group, gens
    using QuantumClifford.ECC: two_block_group_algebra_codes, generalized_bicycle_codes, code_k, code_n, distance, DistanceMIPAlgorithm

    @testset "minimum distance properties: GB" begin
        # [48, 6, 8]] GB code, # minimum distance is exact, d = 8
        l = 24
        c = generalized_bicycle_codes([0, 2, 8, 15], [0, 2, 12, 17], l)
        # minimum distance is exact, d = 8
        i = rand(1:code_k(c))
        # By default, the minimum distance for the Z-type logical operator is computed.
        # The minimum distance for X-type logical operators is the same.
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == distance(c, DistanceMIPAlgorithm(logical_operator_type=:Z; solver=HiGHS)) == 8
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
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == distance(c, DistanceMIPAlgorithm(logical_operator_type=:Z; solver=HiGHS)) == 5
    end
end
