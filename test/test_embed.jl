@testitem "embed PauliOperator" begin
    @test embed(5,3,P"-Y") == P"-__Y__"
    @test embed(5,(3,5),P"-YZ") == P"-__Y_Z"
    @test embed(5,[3,5],P"-YZ") == P"-__Y_Z"
    @test_throws ArgumentError embed(5,[3,5],P"-YZX")
    @test_throws ArgumentError embed(5,3,P"-ZX")
end
