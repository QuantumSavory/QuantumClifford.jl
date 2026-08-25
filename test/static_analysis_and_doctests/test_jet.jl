@testitem "JET checks" tags=[:jet] begin
    using JET
    using Test
    using QuantumClifford
    using QuantumInterface

    rep = JET.report_package(QuantumClifford, target_modules=[QuantumClifford, QuantumInterface])
    @show rep
    @test length(JET.get_reports(rep)) <= 5
    @test_broken length(JET.get_reports(rep)) == 0
end
