@testitem  "ECC GB Codes - slow" tags=[:ecc, :ecc_bespoke_checks] begin
    using Hecke
    using HiGHS
    using JuMP
    using QuantumClifford.ECC.QECCore: code_k, code_n, distance, rate
    using QuantumClifford.ECC: generalized_bicycle_codes, code_k, code_n, DistanceMIPAlgorithm

    @testset "generalized bicycle codes from arXiv:2203.17216" begin
        # codes taken from https://github.com/QEC-pages/GB-codes
        c = generalized_bicycle_codes([0, 2], [0, 1], 5)
        @test code_n(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        c = generalized_bicycle_codes([0, 3], [0, 1], 11)
        @test code_n(c) == 22
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        c = generalized_bicycle_codes([0, 5], [0, 1], 13)
        @test code_n(c) == 26
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5
        c = generalized_bicycle_codes([0, 4], [0, 1], 19)
        @test code_n(c) == 38
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5
        c = generalized_bicycle_codes([0, 8], [0, 1], 29)
        @test code_n(c) == 58
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 7
        c = generalized_bicycle_codes([0, 8], [0, 1], 37)
        @test code_n(c) == 74
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8
        c = generalized_bicycle_codes([0, 8], [0, 1], 53)
        @test code_n(c) == 106
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 9

        c = generalized_bicycle_codes([0, 1, 2, 3], [0, 1], 5)
        @test code_n(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        c = generalized_bicycle_codes([0, 1, 2, 5], [0, 1], 11)
        @test code_n(c) == 22
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5
        c = generalized_bicycle_codes([0, 1, 2, 5], [0, 1], 13)
        @test code_n(c) == 26
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5
        c = generalized_bicycle_codes([0, 2, 3, 7], [0, 1], 19)
        @test code_n(c) == 38
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 7
        c = generalized_bicycle_codes([0, 1, 3, 10], [0, 1], 29)
        @test code_n(c) == 58
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 9
        c = generalized_bicycle_codes([0, 1, 4, 14], [0, 1], 37)
        @test code_n(c) == 74
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 11

        c = generalized_bicycle_codes([0, 1, 2, 3], [0, 1],  5)
        @test code_n(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        c = generalized_bicycle_codes([0, 1, 2, 3], [0, 1],  7)
        @test code_n(c) == 14
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        c = generalized_bicycle_codes([0, 1, 2, 5], [0, 1], 11)
        @test code_n(c) == 22
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5
        c = generalized_bicycle_codes([0, 1, 2, 5], [0, 1], 13)
        @test code_n(c) == 26
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 5
        c = generalized_bicycle_codes([0, 1, 3, 9], [0, 1], 17)
        @test code_n(c) == 34
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 7
        c = generalized_bicycle_codes([0, 1, 2, 8], [0, 1], 19)
        @test code_n(c) == 38
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 7
        c = generalized_bicycle_codes([0, 1, 3, 9], [0, 1], 23)
        @test code_n(c) == 46
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8
        c = generalized_bicycle_codes([0, 1, 3, 10], [0, 1], 29)
        @test code_n(c) == 58
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 9
        c = generalized_bicycle_codes([0, 5, 11,17], [0, 1], 31)
        @test code_n(c) == 62
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 10
        c = generalized_bicycle_codes([0, 1, 4, 14], [0, 1], 37)
        @test code_n(c) == 74
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 11

        c = generalized_bicycle_codes([0, 1, 2, 4, 5, 8], [0, 1], 11)
        @test code_n(c) == 22
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6
        c = generalized_bicycle_codes([0, 2, 3, 6, 7, 9], [0, 1], 13)
        @test code_n(c) == 26
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6
        c = generalized_bicycle_codes([0, 1, 2, 3, 4, 8], [0, 1], 19)
        @test code_n(c) == 38
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 7
    end
end
