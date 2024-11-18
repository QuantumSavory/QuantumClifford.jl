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

##

@testset "inverse sparsity of k-qubit channels" begin
    # In general, the increase in Λ(χ) due to k-qubit channels is limited to a maximum factor of 16ᵏ.
    k = 10
    for n in 1:10
        for repetition in 1:1000
            stab = random_stabilizer(n)
            pauli = random_pauli(n)
            genstab = GeneralizedStabilizer(stab)
            invsparsity_before = genstab |> invsparsity
            i = rand(1:n)
            nc = embed(n, i, pcT)
            for i in 1:k
                apply!(genstab, nc) # in-place
            end
            invsparsity_after = genstab |> invsparsity
            @test (invsparsity_after/invsparsity_before) <= 16^k
        end
    end
end

##

@test_throws ArgumentError GeneralizedStabilizer(S"XX")
@test_throws ArgumentError PauliChannel(((P"X", P"Z"), (P"X", P"ZZ")), (1,2))
@test_throws ArgumentError PauliChannel(((P"X", P"Z"), (P"X", P"Z")), (1,))
@test_throws ArgumentError UnitaryPauliChannel((P"X", P"ZZ"), (1,2))
@test_throws ArgumentError UnitaryPauliChannel((P"X", P"Z"), (1,))
@test_throws MethodError project!(GeneralizedStabilizer(S"X"), P"X")
