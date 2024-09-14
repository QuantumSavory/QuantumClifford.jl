@testitem "JET checks" tags=[:jet] begin
using JET
using ArrayInterface
using Static
using Graphs
using StridedViews
using LinearAlgebra
using AbstractAlgebra, Hecke


rep = report_package("QuantumClifford";
    ignored_modules=(
        AnyFrameModule(Graphs.LinAlg),
        AnyFrameModule(Graphs.SimpleGraphs),
        AnyFrameModule(ArrayInterface),
        AnyFrameModule(Static),
        AnyFrameModule(StridedViews),
        AnyFrameModule(LinearAlgebra),
        AnyFrameModule(Hecke),
        AnyFrameModule(AbstractAlgebra),
    )
)
@show rep
@test_broken length(JET.get_reports(rep)) == 0
@test length(JET.get_reports(rep)) <= 23
end
