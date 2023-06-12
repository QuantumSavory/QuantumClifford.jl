using QuantumClifford
using QuantumClifford: StabMixture, rowdecompose, PauliChannel
using Test
using InteractiveUtils
using Random

##

@testset "Pauli decomposition into destabilizers" begin
    for n in [1,2,63,64,65,66]
        n = 6
        p = random_pauli(n; nophase=true)
        s = random_destabilizer(n)
        phase, b, c = rowdecompose(p,s)
        p0 = zero(p)
        for (i,f) in pairs(b)
            f && mul_left!(p0, destabilizerview(s), i)
        end
        for (i,f) in pairs(c)
            f && mul_left!(p0, stabilizerview(s), i)
        end
        @test (-im)^phase*p0 == p
    end
end

##

@testset "PauliChannel T gate ^4 = Id" begin
    tgate = PauliChannel(
    [(I, I), (I, Z), (Z, I), (Z, Z)],
    [cos(π/8)^2, -im*sin(π/8)*cos(π/8),  im*sin(π/8)*cos(π/8), sin(π/8)^2]
    )

    state = StabMixture(S"-Z")

    apply!(state, tgate)
    apply!(state, tgate)
    apply!(state, tgate)
    apply!(state, tgate)

    @test state.destabweights |> values |> collect == [1]
    @test state.destabweights |> keys |> collect == [([1],[1])]
end
