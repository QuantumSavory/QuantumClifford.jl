@testset "ECC Python decoders" begin
    using Test
    using QuantumClifford
    using QuantumClifford.ECC

    import PyQDecoders
    import PyTesseractDecoder

    include("../test_ecc_base.jl")

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

    @testset "tesseract decoder (tesseract-decoder via PyTesseractDecoder)" begin
        codes = [
            Surface(8, 8),
            Toric(8, 8),
        ]

        noise = 0.02
        setups = [
            CommutationCheckECCSetup(noise),
        ]

        for c in codes
            Hstab = parity_checks(c)
            stab = QuantumClifford.stab_to_gf2(Hstab) # [X | Z]
            s, n2 = size(stab)
            n = n2 ÷ 2
            check_dense = hcat(@view(stab[:, n+1:end]), @view(stab[:, 1:n])) # [Z | X]

            d = TesseractDecoder(c; errorrate=noise, det_beam=50)

            # Basic consistency: the returned error guess matches the input syndrome.
            for col in 1:min(20, 2n)
                syndrome = @view check_dense[:, col]
                guess = QuantumClifford.ECC.decode(d, syndrome)
                predicted_syndrome = falses(s)
                for j in findall(guess)
                    predicted_syndrome .⊻= @view check_dense[:, j]
                end
                @test predicted_syndrome == syndrome
            end

            for setup in setups
                e = evaluate_decoder(d, setup, 2000)
                @test max(e...) <= 0.2/20
            end
        end
    end
end
