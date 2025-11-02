using Test
Base.include(Main, "../src/QuantumClifford.jl")
using .QuantumClifford

@testset "reinterpret" begin
    # PauliOperator roundtrip via reinterpret
    p = P"XYZ"
    p2 = reinterpret(UInt8, p)
    p3 = reinterpret(UInt64, p2)
    @test xbit(p) == xbit(p3)
    @test zbit(p) == zbit(p3)

    # Tableau roundtrip
    t = QuantumClifford.Tableau([P"XX", P"ZZ"])
    t2 = reinterpret(UInt8, t)
    t3 = reinterpret(UInt64, t2)
    for i in 1:length(t.phases)
        @test xbit(t[i]) == xbit(t3[i])
        @test zbit(t[i]) == zbit(t3[i])
    end
end
