@testitem "ECC Decoder" begin
    using QuantumClifford.ECC

    import PyQDecoders
    import LDPCDecoders

    import Nemo: GF
    import LinearAlgebra

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
                    @assert max(e...) < noise/4
                end
            end
        end
    end

    ##

    @testset "belief prop decoders, good for sparse codes" begin
        codes = [
                 # TODO
                ]

        noise = 0.001

        setups = [
                  CommutationCheckECCSetup(noise),
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
                @assert max(e...) < noise/5
            end
        end
    end
end

# TODO add generalized bicycle codes, after which maybe we should remove some of the above codes

other_lifted_product_codes = []

# from https://arxiv.org/abs/2202.01702v3

l = 63
R = PermutationGroupRing(GF(2), l)
A = zeros(R, 7, 7)
x = R(cyclic_permutation(1, l))
A[LinearAlgebra.diagind(A)] .= x^27
A[LinearAlgebra.diagind(A, -1)] .= x^54
A[LinearAlgebra.diagind(A, 6)] .= x^54
A[LinearAlgebra.diagind(A, -2)] .= R(1)
A[LinearAlgebra.diagind(A, 5)] .= R(1)

B = reshape([(1 + x + x^6)'], (1, 1))

push!(other_lifted_product_codes, LPCode(A, B))

@testset "belief prop decoders, good for sparse codes" begin
    codes = vcat(LP04, LP118, other_lifted_product_codes)

    noise = 0.001

    setups = [
        CommutationCheckECCSetup(noise),
        NaiveSyndromeECCSetup(noise, 0),
        ShorSyndromeECCSetup(noise, 0),
    ]
    # lifted product codes currently trigger errors in syndrome circuits

    for c in codes
        for s in setups
            for d in [c -> PyBeliefPropOSDecoder(c, maxiter=10)]
                nsamples = code_n(c) > 400 ? 1000 : 100000
                # take fewer samples for larger codes to save time
                e = evaluate_decoder(d(c), s, nsamples)
                # @show c
                # @show s
                # @show e
                @assert max(e...) < noise / 4 (c, s, e)
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
            @assert max(e...) < noise/5
        end
    end
end
