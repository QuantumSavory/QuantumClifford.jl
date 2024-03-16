using Test
using QuantumClifford
using QuantumClifford.ECC

import PyQDecoders
import LDPCDecoders

@testset "table decoder, good for small codes" begin
    codes = [
        Steane7(),
        Shor9(),
        Perfect5(),
        Cleve8(),
        Gottesman(3),
        Gottesman(4),
        #Gottesman(5), performs poorly
    ]

    noise = 0.001

    setups = [
        CommutationCheckECCSetup(noise/2),
        NaiveSyndromeECCSetup(noise, 0),
        ShorSyndromeECCSetup(noise, 0),
    ]

    for c in codes
        for s in setups
            for d in [TableDecoder]
                e = evaluate_decoder(d(c), s, 100000)
                #@show c
                #@show s
                #@show e
                @assert max(e...) < noise/3
            end
        end
    end
end

##

@testset "belief prop decoders, good for small codes" begin
    codes = [
        #Gottesman(3),
    ]

    noise = 0.001

    setups = [
        CommutationCheckECCSetup(noise/2),
        NaiveSyndromeECCSetup(noise, 0),
        ShorSyndromeECCSetup(noise, 0),
    ]

    for c in codes
        for s in setups
            for d in [c->PyBeliefPropOSDecoder(c, maxiter=10)]
                e = evaluate_decoder(d(c), s, 100000)
                @show c
                @show s
                @show e
                @assert max(e...) < noise/4
            end
        end
    end
end

##

@testset "matching decoder, good as long as column weight of the code is limited" begin
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
