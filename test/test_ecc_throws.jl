@testitem "ECC throws" tags=[:ecc] begin

    using QuantumClifford.ECC: ReedMuller, BCH, RecursiveReedMuller, Golay, Triangular488, Triangular666, Hamming

    @test_throws ArgumentError ReedMuller(-1, 3)
    @test_throws ArgumentError ReedMuller(1, 0)
    @test_throws ArgumentError ReedMuller(4, 2)

    @test_throws ArgumentError BCH(2, 2)
    @test_throws ArgumentError BCH(-2, 2)
    @test_throws ArgumentError BCH(2, -3)
    @test_throws ArgumentError BCH(3, 3)
    @test_throws ArgumentError BCH(3, 4)
    @test_throws ArgumentError BCH(4, 4)
    @test_throws ArgumentError BCH(4, 100)

    @test_throws ArgumentError RecursiveReedMuller(-1, 3)
    @test_throws ArgumentError RecursiveReedMuller(1, 0)
    @test_throws ArgumentError RecursiveReedMuller(4, 2)

    @test_throws ArgumentError Golay(21)

    @test_throws ArgumentError Hamming(1)
  
    @test_throws ArgumentError Triangular488(8)
    @test_throws ArgumentError Triangular488(1)
    @test_throws ArgumentError Triangular666(8)
    @test_throws ArgumentError Triangular666(1)
end
