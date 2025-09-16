@testitem "ECC Syndromes" tags=[:ecc, :ecc_syndrome_measurement] begin
    using QuantumClifford: mul_left!, embed
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")

    using QuantumClifford: Tableau
    reinterpret_frame(frame) = PauliFrame(reinterpret_stab(frame.frame), copy(frame.measurements))
    reinterpret_stab(s) = Stabilizer(Tableau(copy(phases(s)), nqubits(s), collect(reinterpret(UInt8, collect(s.tab.xzs)))[[1:1+(nqubits(s)-1)÷8;end÷2+1:end÷2+1+(nqubits(s)-1)÷8],:]))
    reinterpret_p(p) = PauliOperator(p.phase, nqubits(p), collect(reinterpret(UInt8, p.xz))[[1:1+(nqubits(p)-1)÷8;end÷2+1:end÷2+1+(nqubits(p)-1)÷8]])

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
            pₙ = embed(naive_qubits+1, 1:dataqubits, p) # +1 to account for the buffer qubit hidden in pauli frames
            pₛ = embed(shor_qubits+1, 1:dataqubits, p)  # +1 to account for the buffer qubit hidden in pauli frames
            mul_left!(naive_frames.frame, pₙ)
            mul_left!(shor_frames.frame, pₛ)
            # run the syndrome circuits using the public API
            pftrajectories(naive_frames, naive_scirc)
            pftrajectories(shor_frames, shor_scirc)
            @test pfmeasurements(naive_frames) == pfmeasurements(shor_frames)[:,shor_bits]

            # just for completeness, let's also try bitpacking in UInt8 instead of the default UInt
            _naive_frames = PauliFrame(nframes, naive_qubits, syndromebits)
            _shor_frames = PauliFrame(nframes, shor_qubits, last(shor_bits))
            naive_uint8 = reinterpret_frame(_naive_frames)
            shor_uint8 = reinterpret_frame(_shor_frames)
            pftrajectories(naive_uint8, ecirc)
            pftrajectories(shor_uint8, vcat(ecirc, shor_cat_scirc))
            p_uint8 = reinterpret_p(p)
            pₙ_uint8 = embed(naive_qubits+1, 1:dataqubits, p_uint8)
            pₛ_uint8 = embed(shor_qubits+1, 1:dataqubits, p_uint8)
            mul_left!(naive_uint8.frame, pₙ_uint8)
            mul_left!(shor_uint8.frame, pₛ_uint8)
            pftrajectories(naive_uint8, naive_scirc)
            pftrajectories(shor_uint8, shor_scirc)
            @test pfmeasurements(shor_uint8)[:,shor_bits] == pfmeasurements(shor_frames)[:,shor_bits] == pfmeasurements(naive_frames) == pfmeasurements(naive_uint8)
        end
    end

    @testset "naive and shor measurement circuits" begin
        for (i,c) in enumerate(all_testable_code_instances())
            pframe_naive_vs_shor_syndrome(c)
        end
    end
end
