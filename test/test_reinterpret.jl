using QuantumClifford

@testitem "Reinterpret" begin
    @testset "reinterpret" begin
        # PauliOperator roundtrip via reinterpret
        p = QuantumClifford._P_str("XYZ")
        p2 = QuantumClifford.reinterpret(UInt8, p)
        p3 = QuantumClifford.reinterpret(UInt64, p2)
        @test xbit(p) == xbit(p3)
        @test zbit(p) == zbit(p3)

        # Tableau roundtrip via reinterpret
        t = QuantumClifford.Tableau([QuantumClifford._P_str("XX"), QuantumClifford._P_str("ZZ")])
        t2 = QuantumClifford.reinterpret(UInt8, t)
        t3 = QuantumClifford.reinterpret(UInt64, t2)
        @test QuantumClifford.stab_to_gf2(t) == QuantumClifford.stab_to_gf2(t3)
    end
end
