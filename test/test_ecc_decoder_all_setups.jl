@testitem "ECC Decoder" tags=[:ecc, :ecc_decoding] begin
    using Test
    using QuantumClifford
    using QuantumClifford.ECC

    import PyQDecoders
    import LDPCDecoders
    import Sys

    if !Sys.iswindows()
        import PyTesseractDecoder
    else
        @test_broken false # TODO tesseract-decoder is not available on Windows
    end

    include("test_ecc_base.jl")

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

    @testset "in-place and batch decoder API" begin
        c = Steane7()
        d = CSSTableDecoder(c)
        s = code_s(c)
        n = code_n(c)

        syndrome = falses(s)
        expected = QuantumClifford.ECC.decode(d, syndrome)

        out = trues(2n)
        ret = QuantumClifford.ECC.decode!(out, d, syndrome)
        @test ret === out
        @test out == expected

        syndromes = falses(3, s)
        syndromes[2, 1] = true

        expected_batch = QuantumClifford.ECC.batchdecode(d, syndromes)
        out_batch = trues(3, 2n)
        ret_batch = QuantumClifford.ECC.batchdecode!(out_batch, d, syndromes)
        @test ret_batch === out_batch
        @test out_batch == expected_batch

        out_batch_threaded = fill(true, 3, 2n)
        ret_batch_threaded = QuantumClifford.ECC.batchdecode!(out_batch_threaded, d, syndromes; threaded=true)
        @test ret_batch_threaded === out_batch_threaded
        @test out_batch_threaded == expected_batch

        @test QuantumClifford.ECC.decode_batch(d, syndromes) == expected_batch
        @test QuantumClifford.ECC.decode_batch!(fill(true, 3, 2n), d, syndromes; threaded=true) == expected_batch
    end

@testset "additional edge case coverage for decoders" begin
        c = Steane6()
        s = code_s(c)
        n = code_n(c)
        Hx = parity_matrix_x(c)
        cx = size(Hx, 0)

        td = TableDecoder(c)
        ctd = ClassicalTableDecoder(Matrix{Bool}(Hx))
        csstd = CSSTableDecoder(c)

        syndromes = falses(2, s)
        syndromes[1, 1] = true

        bp_decoder = BeliefPropDecoder(c; maxiter=1)
        out = trues(1n)
        ret = QuantumClifford.ECC.decode!(out, bp_decoder, falses(s))
        @test ret === out

        # Wrong-dimension errors for decode!
        @test_throws ArgumentError QuantumClifford.ECC.decode!(trues(0), td, falses(s))
        @test_throws ArgumentError QuantumClifford.ECC.decode!(trues(0), ctd, falses(cx))
        @test_throws ArgumentError QuantumClifford.ECC.decode!(trues(0), csstd, falses(s))
        @test_throws ArgumentError QuantumClifford.ECC.decode!(trues(1n), csstd, trues(s+1))

        # Wrong-dimension errors for batchdecode! on TableDecoder
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(2,2n), td, trues(3,s+1))
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(1,2n), td, syndromes)
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(2,1), td, syndromes)

        # Wrong-dimension errors for batchdecode! on ClassicalTableDecoder
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(2,7), ctd, trues(3,cx+1))
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(1,7), ctd, syndromes[:, 1:cx])
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(2,1), ctd, syndromes[:, 1:cx])

        # Wrong-dimension errors for batchdecode! on CSSTableDecoder
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(2,2n), csstd, trues(3,s+1))
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(1,2n), csstd, syndromes)
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(2,1), csstd, syndromes)

        @test_throws ArgumentError QuantumClifford.ECC.batchdecode(bp_decoder, trues(2,s+1))
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(2,2n), bp_decoder, trues(3,s+1))
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(1,2n), bp_decoder, syndromes)
        @test_throws ArgumentError QuantumClifford.ECC.batchdecode!(trues(2,1), bp_decoder, syndromes)

        # nothing on unknown syndrome for decode! / decode
        @test isnothing(QuantumClifford.ECC.decode!(trues(1n), td, trues(s)))
        @test isnothing(QuantumClifford.ECC.decode!(trues(n), ctd, trues(cx)))
        @test isnothing(QuantumClifford.ECC.decode(ctd, trues(cx)))

        # CSSTableDecoder.decode! partial failure: X sub-decoder fails first
        x_only_fail = [trues(cx); falses(s-cx)]
        @test isnothing(QuantumClifford.ECC.decode!(trues(1n), csstd, x_only_fail))

        # CSSTableDecoder.decode! partial failure: Z sub-decoder fails (after X succeeds)
        z_only_fail = [falses(cx); trues(s-cx)]
        @test isnothing(QuantumClifford.ECC.decode!(trues(1n), csstd, z_only_fail))

        @test_logs (:warn,) QuantumClifford.ECC.batchdecode!(falses(2,2n), td, syndromes; threaded=true)
        ret = QuantumClifford.ECC.batchdecode!(Matrix{Bool}(undef,2,2n), td, syndromes; threaded=true)
        @test size(ret) == (2, 2n)

        @test_logs (:warn,) QuantumClifford.ECC.batchdecode!(falses(2,n), ctd, syndromes[:, 1:cx]; threaded=true)
        ret = QuantumClifford.ECC.batchdecode!(Matrix{Bool}(undef,2,n), ctd, syndromes[:, 1:cx]; threaded=true)
        @test size(ret) == (2, n)


        @test_logs (:warn,) QuantumClifford.ECC.batchdecode!(falses(2,2n), csstd, syndromes; threaded=true)
        @test_logs (:warn,) QuantumClifford.ECC.batchdecode!(Matrix{Bool}(undef,2,2n), bp_decoder, syndromes; threaded=true)
    end

    if !Sys.iswindows()
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

    end
end
