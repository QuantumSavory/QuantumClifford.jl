using QuantumClifford
using QuantumClifford: GeneralizedStabilizer, rowdecompose, PauliChannel, invsparsity, mul_left!, mul_right!
using Test
using InteractiveUtils
using Random

##

@testset "Pauli decomposition into destabilizers" begin
    for n in [1,2,63,64,65,66,300]
        p = random_pauli(n; nophase=false, realphase=true)
        s = random_destabilizer(n)
        phase, b, c = rowdecompose(p,s)
        p0 = zero(p)
        for (i,f) in pairs(b)
            f && mul_right!(p0, destabilizerview(s), i)
        end
        for (i,f) in pairs(c)
            f && mul_right!(p0, stabilizerview(s), i)
        end
        @test (im)^phase*p0 == p
    end
end

##

@testset "PauliChannel T gate ^4 = Id" begin
    tgate = pcT
    state = GeneralizedStabilizer(S"X")

    apply!(state, tgate)
    apply!(state, tgate)
    apply!(state, tgate)
    apply!(state, tgate)

    @test state.destabweights |> values |> collect ≈ [1]
    @test state.destabweights |> keys |> collect == [([1],[1])]
end

##

@testset "Inverse sparsity" begin
    for n in 1:5
        s = random_stabilizer(n)
        gs = GeneralizedStabilizer(s)
        for i in 1:rand(1:4)
            apply!(gs, embed(n, i, pcT))
        end
        # Λ(χ) ≤ 4ⁿ
        @test invsparsity(gs) <= 4^n
        channels = [embed(n, i, pcT) for i in 1:rand(1:4)]
        # Λ(ϕᵢⱼ) ≤ 4ⁿ
        @test all(invsparsity(channel) <= 4^n for channel in channels)
    end
end

@test_throws ArgumentError GeneralizedStabilizer(S"XX")
@test_throws ArgumentError PauliChannel(((P"X", P"Z"), (P"X", P"ZZ")), (1,2))
@test_throws ArgumentError PauliChannel(((P"X", P"Z"), (P"X", P"Z")), (1,))
@test_throws ArgumentError UnitaryPauliChannel((P"X", P"ZZ"), (1,2))
@test_throws ArgumentError UnitaryPauliChannel((P"X", P"Z"), (1,))
