@testitem "ECC Syndromes" begin
    using QuantumClifford: mul_left!, embed
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")

    function pframe_naive_vs_shor_syndrome(code)
        ecirc = naive_encoding_circuit(code)
        naive_scirc, naive_ancillaries = naive_syndrome_circuit(code)
        shor_cat_scirc, shor_scirc, shor_ancillaries, shor_bits = shor_syndrome_circuit(code)
        nframes = 10
        dataqubits = code_n(code)
        syndromebits = code_s(code)
        naive_qubits = dataqubits + syndromebits
        shor_qubits = dataqubits + shor_ancillaries
        # no noise
        naive_frames = PauliFrame(nframes, naive_qubits, syndromebits)
        shor_frames = PauliFrame(nframes, shor_qubits, last(shor_bits))
        naive_circuit = vcat(ecirc, naive_scirc)
        shor_circuit = vcat(ecirc, shor_cat_scirc, shor_scirc)
        pftrajectories(naive_frames, naive_circuit)
        pftrajectories(shor_frames, shor_circuit)
        @test pfmeasurements(naive_frames) == pfmeasurements(shor_frames)[:,shor_bits]
        # with errors
        for _ in 1:10
            naive_frames = PauliFrame(nframes, naive_qubits, syndromebits)
            shor_frames = PauliFrame(nframes, shor_qubits, last(shor_bits))
            pftrajectories(naive_frames, ecirc)
            pftrajectories(shor_frames, vcat(ecirc, shor_cat_scirc))
            # manually injecting the same type of noise in the frames -- not really a user accessible API
            p = random_pauli(dataqubits, realphase=true)
            pₙ = embed(naive_qubits, 1:dataqubits, p)
            pₛ = embed(shor_qubits, 1:dataqubits, p)
            mul_left!(naive_frames.frame, pₙ)
            mul_left!(shor_frames.frame, pₛ)
            # run the syndrome circuits using the public API
            pftrajectories(naive_frames, naive_scirc)
            pftrajectories(shor_frames, shor_scirc)
            @test pfmeasurements(naive_frames) == pfmeasurements(shor_frames)[:,shor_bits]
        end
    end

    @testset "naive and shor measurement circuits" begin
        for c in all_testablable_code_instances()
            pframe_naive_vs_shor_syndrome(c)
        end
    end
end
