@testitem "ECC GB" begin
    using Hecke
    using HiGHS
    using JuMP
    using QuantumClifford.ECC: generalized_bicycle_codes, code_k, code_n, distance, DistanceMIPAlgorithm

    # codes taken from Table 1 of [lin2024quantum](@cite)
    # Abelian 2BGA codes can be viewed as GB codes.
    @testset "GB codes" begin
        # [[70, 8, 10]]
        c = generalized_bicycle_codes([0, 15, 16, 18], [0, 1, 24, 27], 35)
        @test code_n(c) == 70 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_qubit=1)) == 10
        # [[54, 6, 9]]
        c = generalized_bicycle_codes([0, 1, 3, 7], [0, 1, 12, 19], 27)
        @test code_n(c) == 54 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_qubit=1)) == 9
        # [[60, 6, 10]]
        c = generalized_bicycle_codes([0 , 10, 6, 13], [0, 25, 16, 12], 30)
        @test code_n(c) == 60 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_qubit=1)) == 10
        # [[72, 8, 10]]
        c = generalized_bicycle_codes([0, 9, 28, 31], [0, 1, 21, 34], 36)
        @test code_n(c) == 72 && code_k(c) == 8
        # [[72, 10, 9]]
        c = generalized_bicycle_codes([0, 9, 28, 13], [0, 1, 3, 22], 36)
        @test code_n(c) == 72 && code_k(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_qubit=1)) == 9
    end
end
