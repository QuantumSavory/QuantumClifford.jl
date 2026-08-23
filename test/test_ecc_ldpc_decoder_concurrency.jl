@testitem "LDPC decoders support concurrent reuse" tags=[:ecc, :ecc_decoding] begin
    using Random: Xoshiro, rand
    using QuantumClifford.ECC: BeliefPropDecoder, BitFlipDecoder, Steane7, decode

    syndromes = [rand(Xoshiro(i), Bool, 6) for i in 1:32]
    for decoder_type in (BeliefPropDecoder, BitFlipDecoder)
        decoder = decoder_type(Steane7(); errorrate=0.1, maxiter=10)
        expected = decode.(Ref(decoder), syndromes)
        results = Matrix{Any}(undef, 32, length(syndromes))
        @sync for repetition in axes(results, 1), i in eachindex(syndromes)
            Threads.@spawn results[repetition, i] = decode(decoder, syndromes[i])
        end
        @test all(results[repetition, i] == expected[i]
                  for repetition in axes(results, 1), i in eachindex(syndromes))
    end
end
