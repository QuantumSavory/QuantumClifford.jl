@testitem "Reinterpret" begin
    @testset "reinterpret" begin
        # PauliOperator roundtrip via reinterpret
        p = P"XYZ"
        p2 = reinterpret(UInt8, p)
        p3 = reinterpret(UInt64, p2)
        @test xbit(p) == xbit(p3)
        @test zbit(p) == zbit(p3)

        # Tableau roundtrip via reinterpret
        t = Tableau([P"XX", P"ZZ"])
        t2 = reinterpret(UInt8, t)
        t3 = reinterpret(UInt64, t2)
        @test stab_to_gf2(t) == stab_to_gf2(t3)
    end
end
