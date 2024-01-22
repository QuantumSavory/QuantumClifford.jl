using Test
using QuantumClifford
using QuantumClifford.ECC
import PyQDecoders


@testset "table decoder" begin
    codes = [
        Steane7(),
        Shor9(),
        Perfect5(),
        Cleve8(),
    ]

    noise = 0.001

    setups = [
        CommutationCheckECCSetup(noise/2),
        NaiveSyndromeECCSetup(noise, 0),
        ShorSyndromeECCSetup(noise, 0),
    ]

    for c in codes
        for s in setups
            e = evaluate_decoder(TableDecoder(c), s, 100000)
            #@show c
            #@show s
            #@show e
            @assert max(e...) < noise/4
        end
    end
end

##

@testset "matching decoder" begin
    codes = [
        Toric(8,8),
        Toric(9,9)
    ]

    noise = 0.01

    setups = [
        CommutationCheckECCSetup(noise/2),
        NaiveSyndromeECCSetup(noise, 0),
        ShorSyndromeECCSetup(noise, 0),
    ]

    for c in codes
        for s in setups
            e = evaluate_decoder(PyMatchingDecoder(c), s, 10000)
            #@show c
            #@show s
            #@show e
            @assert max(e...) < noise/5
        end
    end
end
