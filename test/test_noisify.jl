@testitem "noisify" begin
    using QuantumClifford

    @testset "simple API" begin
        result = noisify([sHadamard(1)], PauliNoise(1e-3, 1e-3, 1e-3))
        @test length(result) == 2
        @test result[2] == sHadamard(1)
    end

    @testset "errors when idle requested without nqubits" begin
        @test_throws ErrorException noisify([sHadamard(1)], DefaultNoiseModel())
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
        result = noisify([op], PauliNoise(1e-3, 1e-3, 1e-3))
        @test length(result) == 1
        @test result[1] === op
    end
end