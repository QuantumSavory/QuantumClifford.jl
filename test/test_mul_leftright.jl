using QuantumClifford
using QuantumClifford: mul_left!, mul_right!

using Test

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

@testset "Inner product between stabilizer states" begin
    for n in test_sizes
        for _ in 1:20
            p1 = random_pauli(n)
            p2 = random_pauli(n)
            s = random_stabilizer(n)
            i = rand(1:n)
            @test p1*p2 == mul_left!(copy(p2), p1)
            @test p1*p2 == mul_right!(copy(p1), p2)
            @test mul_left!(copy(p2), p1) == (-1)^comm(p1,p2) * mul_right!(copy(p2), p1)
            @test mul_left!(copy(p2), s[i]) == mul_left!(copy(p2), s, i) == s[i]*p2
            @test mul_right!(copy(p2), s[i]) == mul_right!(copy(p2), s, i) == p2*s[i]
            @test mul_left!(copy(s), p2)[i] == p2*s[i]
            @test mul_right!(copy(s), p2)[i] == s[i]*p2
        end
    end
end
