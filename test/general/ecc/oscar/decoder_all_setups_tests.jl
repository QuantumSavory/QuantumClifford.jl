@testset "ECC Decoder" begin
    using Test
    using QuantumClifford
    using QuantumClifford.ECC

    import LDPCDecoders

    include("../test_ecc_base.jl")

    @testset "table decoder, good for small codes" begin
        codes = [
                 all_testable_code_instances(;maxn=10)...
                ]

        noise = 0.001

        setups = [
                  CommutationCheckECCSetup(noise),
                  NaiveSyndromeECCSetup(noise, 0),
                  ShorSyndromeECCSetup(noise, 0),
                 ]

        for c in codes
            for s in setups
                for d in [TableDecoder, CSSTableDecoder]
                    d === CSSTableDecoder && !iscss(c) && continue
                    e = evaluate_decoder(d(c), s, 100000)
                    #@show c
                    #@show s
                    #@show e
                    @test max(e...) < noise/4
                end
            end
        end
    end

    @testset "BitFlipDecoder decoder, good for sparse codes" begin
        codes = [
                 QuantumReedMuller(3),
                 QuantumReedMuller(4)
                ]

        noise = 0.001

        setups = [
                  CommutationCheckECCSetup(noise),
                  NaiveSyndromeECCSetup(noise, 0),
                  ShorSyndromeECCSetup(noise, 0),
                 ]

        for c in codes
            for s in setups
                for d in [c->BitFlipDecoder(c, maxiter=10)]
                    e = evaluate_decoder(d(c), s, 100000)
                    #@show c
                    #@show s
                    #@show e
                    @test max(e...) < noise/4
                end
            end
        end
    end

end
