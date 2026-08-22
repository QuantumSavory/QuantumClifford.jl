@testitem "Tile2D public API" tags=[:ecc, :ecc_base] begin
    using Test
    using QuantumClifford
    using QuantumClifford.ECC

    tile = Tile2D(2, [(0,0)], [(0,0)], 1, 2)
    @test :Tile2D in names(QuantumClifford.ECC)
    @test code_n(tile) == 4
end
