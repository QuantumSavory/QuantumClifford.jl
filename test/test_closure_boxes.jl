@testitem "Closure boxes" tags=[:julia_1_14_required] begin
    using Hecke
    hecke_ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    @test !isnothing(hecke_ext)
    @test isempty(Test.detect_closure_boxes(QuantumClifford, hecke_ext))
end
