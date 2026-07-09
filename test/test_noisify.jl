@testitem "noisify" begin
    using QuantumClifford:AbstractNoiseOp, Reset

    @testset "simple noise insertion" begin
        reset_state = one(Stabilizer, 1)

        circuit = Any[
            sHadamard(1),
            sCNOT(1, 2),
            sMZ(1, 1),
            Reset(reset_state, [2]),
        ]

        noisy = noisify(circuit, PauliNoise(1e-3, 1e-3, 1e-3))

        @test length(noisy) == 8

        @test noisy[1] isa AbstractNoiseOp
        @test noisy[2] == circuit[1]

        @test noisy[3] isa AbstractNoiseOp
        @test noisy[4] == circuit[2]

        @test noisy[5] isa AbstractNoiseOp
        @test noisy[6] == circuit[3]

        @test noisy[7] isa Reset
        @test noisy[8] isa AbstractNoiseOp

        original = copy(circuit)

        noisy = noisify(circuit, PauliNoise(1e-3, 1e-3, 1e-3))

        @test length(circuit) == length(original)
        @test all(circuit[i] === original[i] for i in eachindex(circuit))
    end

    @testset "Empty CircuitNoise leaves circuit unchanged" begin
        circuit = [
            sHadamard(1),
            sCNOT(1, 2),
            sMZ(1, 1),
        ]

        @test noisify(circuit, CircuitNoise()) == circuit
    end

    @testset "Certain ops pass through unchanged" begin
        verify = VerifyOp(one(Stabilizer, 1), [1])
        xor = ClassicalXOR(1, 1)
        existing_noise = NoiseOp(PauliNoise(0.1, 0.1, 0.1), [1])

        circuit = Any[
            verify,
            xor,
            existing_noise,
        ]

        noisy = noisify(circuit, PauliNoise(1e-3, 1e-3, 1e-3))

        @test noisy == circuit
    end
    @testset "skipped ops do not trigger idle noise" begin
        verify = VerifyOp(one(Stabilizer, 1), [1])

        circuit = [
            sHadamard(1),
            verify,
        ]

        model = CircuitNoise(
            idle_noise = PauliNoise(1e-5, 1e-5, 1e-5),
        )

        noisy = noisify(circuit, model)

        @test noisy == circuit
    end
    @testset "idle noise insertion with parallel operations" begin
        circuit = [
            sHadamard(1),
            sHadamard(2),
            sCNOT(1, 3),
            sHadamard(1),
            sMZ(1, 1),
        ]

        model = CircuitNoise(
            idle_noise = PauliNoise(1e-5, 1e-5, 1e-5),
        )

        noisy = noisify(circuit, model)

        filtered = filter(op -> !(op isa AbstractNoiseOp), noisy)
        @test filtered == circuit

        noise_ops = filter(op -> op isa AbstractNoiseOp, noisy)
        @test !isempty(noise_ops)

        q1_ops = filter(op -> Tuple(affectedqubits(op)) == (1,), noise_ops)
        @test isempty(q1_ops)

        q2_ops = filter(op -> Tuple(affectedqubits(op)) == (2,), noise_ops)
        @test length(q2_ops) == 3

        q3_ops = filter(op -> Tuple(affectedqubits(op)) == (3,), noise_ops)
        @test length(q3_ops) == 3

        @test length(noise_ops) == 6
        @test length(noisy) == length(circuit) + length(noise_ops)
    end

    @testset "original circuit order remains" begin

        circuit = Any[
            sHadamard(1),
            sCNOT(1, 2),
            VerifyOp(one(Stabilizer, 1), [1]),
            sMZ(1, 1),
        ]

        model = CircuitNoise(
            single_qubit = PauliNoise(1e-4, 1e-4, 1e-4),
            two_qubit    = PauliNoise(1e-3, 1e-3, 1e-3),
            measurement  = PauliNoise(2e-3, 2e-3, 2e-3),
            idle_noise = PauliNoise(1e-5, 1e-5, 1e-5)
        )

        noisy = noisify(circuit, model)

        filtered = filter(op -> !(op isa AbstractNoiseOp), noisy)

        @test filtered == circuit
    end
end
