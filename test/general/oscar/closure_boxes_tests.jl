@testset "Oscar extension closure boxes" begin
    @static if VERSION >= v"1.14"
        using Oscar
        oscar_ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
        @test !isnothing(oscar_ext)
        @test isempty(Test.detect_closure_boxes(oscar_ext))
    end
end
