@testitem "Inner product between stabilizer states" begin
    using LinearAlgebra
    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    for n in rand(test_sizes)
        s1 = one(Stabilizer, n, basis=:X)
        s2 = one(MixedDestabilizer, n, n)
        @test dot(s1,s2)≈2.0^(-n/2)
        c = random_clifford(n)
        s3 = random_stabilizer(n)
        s4 = random_destabilizer(n)
        @test logdot(s3,s4)==logdot(c*s3,c*s4)
    end
    sa = S" XX
    ZZ"
    sb = S" XZ
    ZX"
    sc = S" XX
    -ZZ"
    @test isnothing(logdot(sa,sb))
    @test isnothing(logdot(S"Y",S"-Y"))
    @test logdot(S"Z",S"X") == 1
    @test logdot(S"X",S"X") == 0
    @test isnothing(logdot(S"X"⊗sa,S"X"⊗sb))
    @test logdot(S"X"⊗sb,S"X"⊗sc) == 1
    @test logdot(S"X"⊗sb,S"Z"⊗sc) == 2
end
