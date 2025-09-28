@testitem "ECC - naive syndrome circuits" tags=[:ecc, :ecc_syndrome_measurement_correctness] begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")

    codes = all_testable_code_instances(; maxn=100)

    function test_naive_syndrome(c::AbstractECC, e::Bool)
        # create a random logical state
        unencoded_qubits = random_stabilizer(code_k(c))
        bufferqubits = one(Stabilizer, code_n(c) - code_k(c))
        logicalqubits = bufferqubits⊗unencoded_qubits
        mctrajectory!(logicalqubits, naive_encoding_circuit(c))
        if e
            #add some noise to logicalqubits
            apply!(logicalqubits, P"X", rand(1:code_n(c)))
            apply!(logicalqubits, P"Z", rand(1:code_n(c)))
        end
        # measure using `project!`
        s1 = copy(logicalqubits)
        syndrome1 = [project!(s1, check)[3] for check in parity_checks(c)]
        # measure using `naive_syndrome_circuit`
        naive_circuit, _ = naive_syndrome_circuit(c)
        ancillaryqubits = one(Stabilizer,code_s(c))
        s2 = copy(logicalqubits)
        syndrome2 = Register(s2⊗ancillaryqubits, falses(code_s(c)))
        mctrajectory!(syndrome2, naive_circuit)
        if !e
            @test all(syndrome1 .== 0)
            @test all(bitview(syndrome2) .== 0)
        end
        @test bitview(syndrome2) == syndrome1.÷2
    end

    @testset "naive syndrome circuits - zero syndrome for logical states" begin
        for c in codes, _ in 1:10
            test_naive_syndrome(c, false)
            test_naive_syndrome(c, true)
        end
    end

    ##

    function test_with_pframes(code)
        ecirc = naive_encoding_circuit(code)
        scirc, _ = naive_syndrome_circuit(code)
        nframes = 10
        dataqubits = code_n(code)
        ancqubits = code_s(code)
        regbits = ancqubits
        frames = PauliFrame(nframes, dataqubits+ancqubits, regbits)
        circuit = vcat(ecirc, scirc)
        pftrajectories(frames, circuit)
        @test sum(pfmeasurements(frames)) == 0
    end

    @testset "naive syndrome circuits - zero syndrome for logical states" begin
        for c in codes, _ in 1:2
            test_with_pframes(c)
        end
    end
end
