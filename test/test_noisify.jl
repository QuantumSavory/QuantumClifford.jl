@testitem "noisify" begin
    using QuantumClifford

    @testset "simple API" begin
        result = noisify([sHadamard(1)], PauliNoise(1e-3, 1e-3, 1e-3))
        @test length(result) == 2
        @test result[2] == sHadamard(1)
    end

    @testset "errors when idle requested without nqubits" begin
        @test_throws "nqubits must be provided" noisify(
            [sHadamard(1)], DefaultNoiseModel()
        )
    end

    @testset "noise inserted before op with correct qubits" begin
        result = noisify([sCNOT(1, 2)], PauliNoise(1e-3, 1e-3, 1e-3))
        @test length(result) == 2
        @test result[1] isa NoiseOp
        @test result[2] == sCNOT(1, 2)
        @test affectedqubits(result[1]) == (1, 2)
    end

    @testset "unknown operations pass through unchanged" begin
        op = ClassicalXOR((1, 2), 3)
        result = noisify([op], PauliNoise(1e-3, 1e-3, 1e-3); nqubits=3)
        @test length(result) == 1
        @test result[1] === op
    end
    @testset "noisify dispatches NoiseModel fields correctly" begin
        model = NoiseModel(
            single_qubit = PauliNoise(0.01, 0.0,  0.0),
            two_qubit    = PauliNoise(0.0,  0.02, 0.0),
            measurement  = PauliNoise(0.0,  0.0,  0.03),
            reset        = PauliNoise(0.04, 0.04, 0.04),
            idle         = PauliNoise(0.05, 0.0,  0.0),
        )

        circuit = [
            sHadamard(1),
            sCNOT(1, 2),
            sX(2),
            sMZ(1),
            sMRZ(3),
        ]

        noisified = noisify(circuit, model; nqubits=3)
        @test length(noisified) == 15

        # sHadamard(1): idle on q2,q3; single-qubit noise on q1
        @testset "sHadamard(1)" begin
            @test noisified[1].noise == model.idle
            @test noisified[1].indices == (2, 3)
            @test noisified[2].noise == model.single_qubit
            @test noisified[2].indices == (1,)
            @test noisified[3] == circuit[1]
        end

        # sCNOT(1,2): idle on q3; two-qubit noise on q1,q2
        @testset "sCNOT(1,2)" begin
            @test noisified[4].noise == model.idle
            @test noisified[4].indices == (3,)
            @test noisified[5].noise == model.two_qubit
            @test noisified[5].indices == (1, 2)
            @test noisified[6] == circuit[2]
        end

        # sX(2): idle on q1,q3; single-qubit noise on q2
        @testset "sX(2)" begin
            @test noisified[7].noise == model.idle
            @test noisified[7].indices == (1, 3)
            @test noisified[8].noise == model.single_qubit
            @test noisified[8].indices == (2,)
            @test noisified[9] == circuit[3]
        end

        # sMZ(1): idle on q2,q3; measurement noise on q1
        @testset "sMZ(1)" begin
            @test noisified[10].noise == model.idle
            @test noisified[10].indices == (2, 3)
            @test noisified[11].noise == model.measurement
            @test noisified[11].indices == (1,)
            @test noisified[12] == circuit[4]
        end

        # sMRZ(3): idle on q1,q2; reset noise on q3
        @testset "sMRZ(3)" begin
            @test noisified[13].noise == model.idle
            @test noisified[13].indices == (1, 2)
            @test noisified[14].noise == model.reset
            @test noisified[14].indices == (3,)
            @test noisified[15] == circuit[5]
        end
    end
end