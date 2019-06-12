using SimpleClifford, Test

function tests()
@testset "Parsing and constructors" begin
    @test P"-iXYZ" == PauliOperator(0x3, 3, vcat(BitArray([1,1,0]).chunks, BitArray([0,1,1]).chunks))
    @test P"-iXYZ" == PauliOperator(0x3, Bool[1,1,0], Bool[0,1,1])
end

@testset "Elementary operations" begin
    @test P"X"*P"Z" == P"-iY"
end

@testset "Stabilizer canonicalization" begin
    s = S"""- XZZZZ_____
            - _YZY___YX_
            - __XXZ__YX_
            + Z_Z_Y__YXZ
            + _____Z____
            + __________
            + __________
            + ______YYX_
            + __________
            + __________"""
    canonicalize!(s)
    t = S"""- XZZZZ_____
            - _YZY___YX_
            - __XXZ__YX_
            + Z_Z_Y__YXZ
            + ______YYX_
            + _____Z____
            + __________
            + __________
            + __________
            + __________"""
    @test s == t
end

end

tests()
