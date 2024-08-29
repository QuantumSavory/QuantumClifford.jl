@testitem "ECC throws" begin
    using QuantumClifford.ECC: ReedMuller, BCH

    @test_throws ArgumentError ReedMuller(-1, 3)
    @test_throws ArgumentError ReedMuller(1, 0)

    @test_throws ArgumentError BCH(2, 2)
    @test_throws ArgumentError BCH(3, 3)
    @test_throws ArgumentError BCH(3, 4)
end
