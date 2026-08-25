@testitem "Oscar extension closure boxes" tags=[:oscar_required, :julia_1_14_required] begin
    using Oscar
    oscar_ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    @test !isnothing(oscar_ext)
    @test isempty(Test.detect_closure_boxes(oscar_ext))
end
