@testitem "noisify" begin
    using QuantumClifford
    using QuantumClifford: AbstractNoiseOp, Reset

    @testset "simple noise model" begin
        s = one(Stabilizer, 2)

        circuit = Any[
            sHadamard(1),
            sCNOT(1, 2),
            sMZ(1, 1),
            Reset(s, [2]),
        ]

        noise = PauliNoise(1e-3, 1e-3, 1e-3)
        noisy = noisify(circuit, noise)

        @test length(noisy) == 8

        @test noisy[1] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[1])) == (1,)
        @test noisy[2] == circuit[1]

        @test noisy[3] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[3])) == (1, 2)
        @test noisy[4] == circuit[2]

        @test noisy[5] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[5])) == (1,)
        @test noisy[6] == circuit[3]


        @test noisy[7] == circuit[4]
        @test noisy[8] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[8])) == (2,)

        @test length(circuit) == 4
        @test circuit[1] == sHadamard(1)
        @test circuit[2] == sCNOT(1, 2)
        @test circuit[3] == sMZ(1, 1)
        @test circuit[4] isa Reset
    end

    @testset "NoNoise leaves circuit unchanged" begin
        circuit = [
            sHadamard(1),
            sCNOT(1, 2),
            sMZ(1, 1),
        ]

        noisy = noisify(circuit, NoNoise())

        @test noisy == circuit
    end

    @testset "CircuitNoise inserts location-specific noise" begin
        circuit = [
            sHadamard(1),
            sHadamard(2),
            sCNOT(1, 3),
            sMZ(1, 1),
        ]

        noise_model = CircuitNoise(
            single_qubit = PauliNoise(1e-4, 1e-4, 1e-4),
            two_qubit    = PauliNoise(1e-3, 1e-3, 1e-3),
            idle_noise   = PauliNoise(1e-5, 1e-5, 1e-5),
            measurement  = PauliNoise(2e-3, 2e-3, 2e-3),
        )

        noisy = noisify(circuit, noise_model)

        @test noisy[1] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[1])) == (3,)

        @test noisy[2] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[2])) == (1,)
        @test noisy[3] == circuit[1]

        @test noisy[4] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[4])) == (2,)
        @test noisy[5] == circuit[2]

        @test noisy[6] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[6])) == (2,)

        @test noisy[7] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[7])) == (1, 3)
        @test noisy[8] == circuit[3]

        @test noisy[end-1] isa AbstractNoiseOp
        @test Tuple(affectedqubits(noisy[end-1])) == (1,)
        @test noisy[end] == circuit[4]
    end
end
