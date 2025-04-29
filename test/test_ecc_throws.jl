@testitem "ECC throws" begin

    using QuantumClifford.ECC: ReedMuller, BCH, RecursiveReedMuller, Golay

    @test_throws ArgumentError ReedMuller(-1, 3)
    @test_throws ArgumentError ReedMuller(1, 0)
    @test_throws ArgumentError ReedMuller(4, 2)

    @test_throws ArgumentError BCH(2, 2)
    @test_throws ArgumentError BCH(3, 4)

    @test_throws ArgumentError RecursiveReedMuller(-1, 3)
    @test_throws ArgumentError RecursiveReedMuller(1, 0)
    @test_throws ArgumentError RecursiveReedMuller(4, 2)

    @test_throws ArgumentError Golay(21)
end
