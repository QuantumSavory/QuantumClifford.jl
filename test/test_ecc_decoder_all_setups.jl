using Test
using QuantumClifford
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

# test codes LP04 and LP118 are from https://arxiv.org/pdf/2111.07029

B04 = Dict(
    7 => [ 0 0 0 0; 0 1 2 5; 0 6 3 1],
    9 => [0 0 0 0; 0 1 6 7; 0 4 5 2],
    17 => [0 0 0 0; 0 1 2 11; 0 8 12 13],
    19 => [0 0 0 0; 0 2 6 9; 0 16 7 11]
)

B118 = Dict(
    16 => [0 0 0 0 0; 0 2 4 7 11; 0 3 10 14 15],
    21 => [0 0 0 0 0; 0 4 5 7 17; 0 14 18 12 11],
    30 => [0 0 0 0 0; 0 2 14 24 25; 0 16 11 14 13],
)

LP04 = []
for l in keys(B04)
    A = map(B04[l]) do x
        (PermutationGroupRing(GF(2), l))(cyclic_permutation(x, l))
    end
    push!(LP04, LPCode(LiftedCode(A), LiftedCode(A)))
end

LP118 = []
for l in keys(B118)
    A = map(B118[l]) do x
        (PermutationGroupRing(GF(2), l))(cyclic_permutation(x, l))
    end
    push!(LP118, LPCode(LiftedCode(A), LiftedCode(A)))
end

other_lifted_product_codes = []

# from https://arxiv.org/abs/2202.01702v3

l = 63
R = PermutationGroupRing(GF(2), l)
A = zeros(R, 7,7)
A[LinearAlgebra.diagind(A)] .= R(cyclic_permutation(27, l))
A[LinearAlgebra.diagind(A, -1)] .= R(cyclic_permutation(54, l))
A[LinearAlgebra.diagind(A, 6)] .= R(cyclic_permutation(54, l))
A[LinearAlgebra.diagind(A, -2)] .= R(1)
A[LinearAlgebra.diagind(A, 5)] .= R(1)

B = zeros(R, 1, 1)
B[1,1] = (R(1) + R(cyclic_permutation(1, l)) + R(cyclic_permutation(6, l)))'

push!(other_lifted_product_codes, LPCode(LiftedCode(A), LiftedCode(B)))


@testset "belief prop decoders, good for sparse codes" begin
    codes = vcat(LP04, LP118, other_lifted_product_codes)

    noise = 0.001

    setups = [
        CommutationCheckECCSetup(noise),
        # NaiveSyndromeECCSetup(noise, 0),
        # ShorSyndromeECCSetup(noise, 0),
    ]
    # lifted product codes currently trigger errors in syndrome circuits

    for c in codes
        for s in setups
            for d in [c->PyBeliefPropOSDecoder(c, maxiter=10)]
                nsamples = code_n(c)>400 ? 1000 : 100000
                # take fewer samples for larger codes to save time
                e = evaluate_decoder(d(c), s, nsamples)
                # @show c
                # @show s
                # @show e
                @assert max(e...) < noise/4 (c, s, e)
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
