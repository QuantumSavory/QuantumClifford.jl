using QuantumClifford

@testitem "Reinterpret" begin
    p = P"XYZ"
    p2 = reinterpret(UInt8, p)
    p3 = reinterpret(UInt64, p2)
    @test xbit(p) == xbit(p3)
    @test zbit(p) == zbit(p3)

    # Tableau roundtrip via reinterpret
    t = T"XX ZZ"
    t2 = reinterpret(UInt8, t)
    t3 = reinterpret(UInt64, t2)
    @test QuantumClifford.stab_to_gf2(t) == QuantumClifford.stab_to_gf2(t3)
end

@testitem "reinterpret edge cases" begin
    @test_throws ArgumentError reinterpret(UInt128, P"XXXXX")

    small_xz = UInt8[0x0, 0x0, 0x0, 0x0]
    p_small = QuantumClifford.PauliOperator(0x0, 2, small_xz)
    @test_throws ArgumentError reinterpret(UInt32, p_small)
end

@testset "reinterpret combinations" begin
    unsigned_types = (UInt8, UInt16, UInt32, UInt64, UInt128)
    ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]

    for n in ns
        for Ti in unsigned_types, Tf in unsigned_types
            len = QuantumClifford._nchunks(n, Ti)
            xz = rand(Ti, len)
            p = PauliOperator(0x0, n, xz)
            try
                p2 = reinterpret(Tf, p)
                p3 = reinterpret(Ti, p2)
                @test xbit(p) == xbit(p3)
                @test zbit(p) == zbit(p3)
            catch e
                @test isa(e, ArgumentError)
            end
        end
    end
end
