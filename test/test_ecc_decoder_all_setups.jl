@testitem "ECC Decoder" tags=[:ecc] begin
    using QuantumClifford.ECC

    import PyQDecoders
    import LDPCDecoders

    include("test_ecc_base.jl")

    @testset "table decoder, good for small codes" begin
        codes = [
                 all_testablable_code_instances(;maxn=10)...
                ]

        noise = 0.001

        setups = [
                  CommutationCheckECCSetup(noise),
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
                    @test max(e...) < noise/4
                end
            end
        end
    end

    ##

    @testset "belief prop decoders, good for sparse codes" begin
        codes = vcat(LP04, LP118, test_gb_codes, other_lifted_product_codes)

        noise = 0.001

        setups = [
            CommutationCheckECCSetup(noise),
            NaiveSyndromeECCSetup(noise, 0),
            ShorSyndromeECCSetup(noise, 0),
        ]

        for c in codes
            for s in setups
                for d in [c -> PyBeliefPropOSDecoder(c, maxiter=2)]
                    nsamples = 10000
                    if true
                        @test_broken false # TODO these are too slow to test in CI
                        continue
                    end
                    e = evaluate_decoder(d(c), s, nsamples)
                    # @show c
                    # @show s
                    # @show e
                    @test max(e...) <= noise
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

    ##

    using Test
    using QuantumClifford
    using QuantumClifford.ECC

    import PyQDecoders
    import LDPCDecoders

    @testset "matching decoder, good as long as column weight of the code is limited" begin
        codes = [
                 Toric(8,8),
                 Toric(9,9),
                 Surface(8,8),
                 Surface(9,9)
                ]

        noise = 0.01

        setups = [
                  CommutationCheckECCSetup(noise),
                  NaiveSyndromeECCSetup(noise, 0),
                  ShorSyndromeECCSetup(noise, 0),
                 ]

        for c in codes
            for s in setups
                e = evaluate_decoder(PyMatchingDecoder(c), s, 10000)
                #@show c
                #@show s
                #@show e
                @test max(e...) < noise/5
            end
        end
    end
end
