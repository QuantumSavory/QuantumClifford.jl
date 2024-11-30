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
    using JuMP
    using GLPK

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
            AnyFrameModule(JuMP),
            AnyFrameModule(GLPK),
    ))

    @show rep
    @show length(JET.get_reports(rep))
    @test length(JET.get_reports(rep)) == 0
end
