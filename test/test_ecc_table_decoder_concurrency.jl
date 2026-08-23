@testitem "Table decoders support concurrent reuse" tags=[:ecc, :ecc_decoding] begin
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: ClassicalTableDecoder, decode

    decoder = TableDecoder(Steane7())
    checks = parity_checks(Steane7())
    syndromes = [comm(single_x(7, qubit), checks) for qubit in 1:7]
    expected = decode.(Ref(decoder), syndromes)
    results = Matrix{Any}(undef, 32, length(syndromes))
    @sync for repetition in axes(results, 1), i in eachindex(syndromes)
        Threads.@spawn results[repetition, i] = decode(decoder, syndromes[i])
    end
    @test all(results[repetition, i] == expected[i]
              for repetition in axes(results, 1), i in eachindex(syndromes))

    classical = ClassicalTableDecoder(Bool[1 1 0; 0 1 1])
    classical_syndromes = (Bool[true, false], Bool[false, true], Bool[true, true])
    classical_expected = decode.(Ref(classical), classical_syndromes)
    classical_results = Matrix{Any}(undef, 32, length(classical_syndromes))
    @sync for repetition in axes(classical_results, 1), i in eachindex(classical_syndromes)
        Threads.@spawn classical_results[repetition, i] = decode(classical, classical_syndromes[i])
    end
    @test all(classical_results[repetition, i] == classical_expected[i]
              for repetition in axes(classical_results, 1), i in eachindex(classical_syndromes))
end
