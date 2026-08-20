@testitem "Aqua" tags=[:aqua] begin
    using Aqua
    Aqua.test_all(QuantumClifford)
    if isdefined(Test, :detect_closure_boxes)
        @test isempty(Test.detect_closure_boxes(QuantumClifford))
    end
end
