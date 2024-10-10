@testitem "JET checks" tags=[:jet] begin
    using JET
    using Test
    using ArrayInterface
    using Static
    using Graphs
    using StridedViews
    using LinearAlgebra
    using Nemo
    using AbstractAlgebra
    using Hecke
    using StaticArrays

    rep = report_package("QuantumClifford";
        ignored_modules=(
            AnyFrameModule(Graphs.LinAlg),
            AnyFrameModule(Graphs.SimpleGraphs),
            AnyFrameModule(ArrayInterface),
            AnyFrameModule(Static),
            AnyFrameModule(StridedViews),
            AnyFrameModule(LinearAlgebra),
            AnyFrameModule(Nemo),
            AnyFrameModule(AbstractAlgebra),
            AnyFrameModule(Hecke),
            AnyFrameModule(StaticArrays),
    ))

    @show rep
    @show length(JET.get_reports(rep))
    @test_broken length(JET.get_reports(rep)) == 0
    @test length(JET.get_reports(rep)) <= 1
end
