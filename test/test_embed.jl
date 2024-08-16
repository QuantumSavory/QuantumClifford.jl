@testitem "embed methods" begin

    @testset "embed PauliOperator" begin
        @test embed(5,3,P"-Y") == P"-__Y__"
        @test embed(5,(3,5),P"-YZ") == P"-__Y_Z"
        @test embed(5,[3,5],P"-YZ") == P"-__Y_Z"
        @test_throws ArgumentError embed(5,[3,5],P"-YZX")
        @test_throws ArgumentError embed(5,3,P"-ZX")
    end

    @testset "embed Stabilizer" begin
        @test embed(5,3,S"-Y") == S"-__Y__"
        @test embed(5,(3,5),S"-YZ") == S"-__Y_Z"
        @test embed(5,[3,5],S"-YZ") == S"-__Y_Z"
        @test embed(5, 5, S"-Y Z") == S"-____Y ____Z"
        @test embed(5, (2,4), S"XX -YZ") == S"_X_X_ -_Y_Z_"
        @test_throws ArgumentError embed(5,[3,5],S"-YZX")
        @test_throws ArgumentError embed(5,3,S"-ZX")
    end

    @testset "embed Tableau" begin
        @test embed(5,3,T"-Y") == T"-__Y__"
        @test embed(5,(3,5),T"-YZ") == T"-__Y_Z"
        @test embed(5,[3,5],T"-YZ") == T"-__Y_Z"
        @test embed(5, 5, T"-Y Z") == T"-____Y ____Z"
        @test embed(5, (2,4), T"XX -YZ") == T"_X_X_ -_Y_Z_"
        @test_throws ArgumentError embed(5,[3,5],T"-YZX")
        @test_throws ArgumentError embed(5,3,T"-ZX")
    end
end
