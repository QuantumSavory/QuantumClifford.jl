@testitem "ECC Color Codes" begin
    using QECCore
    using Test
    using Hecke
    using JuMP
    using HiGHS
    using Hecke: group_algebra, GF, abelian_group, gens
    import QuantumClifford.ECC: distance, DistanceMIPAlgorithm

    @testset "Minimum Distance test of Color Codes" begin
        d = [3, 5, 7, 9, 11]
        code_types = [Triangular488, Triangular666]
        for i in d
            for CodeType in code_types
                code_instance = CodeType(i)
                @test distance(code_instance, DistanceMIPAlgorithm(solver=HiGHS)) == i
            end
        end
    end
end
