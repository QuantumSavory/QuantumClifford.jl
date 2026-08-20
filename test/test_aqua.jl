@testitem "Aqua" tags=[:aqua] begin
    using Aqua
    Aqua.test_all(QuantumClifford)
end

@testitem "Closure boxes" tags=[:aqua] begin
    using Hecke
    if isdefined(Test, :detect_closure_boxes)
        hecke_ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
        @test !isnothing(hecke_ext)
        @test isempty(Test.detect_closure_boxes(QuantumClifford, hecke_ext))
    end
end

@testitem "Oscar extension closure boxes" tags=[:aqua, :oscar_required] begin
    using Oscar
    if isdefined(Test, :detect_closure_boxes)
        oscar_ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
        @test !isnothing(oscar_ext)
        @test isempty(Test.detect_closure_boxes(oscar_ext))
    end
end
