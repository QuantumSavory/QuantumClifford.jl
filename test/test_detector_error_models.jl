@testitem "Stim detector error model import" begin
    using Random

    @testset "basic detector and logical sampling" begin
        model = detector_error_model_circuit(IOBuffer("""
        detector D0
        detector D1
        logical_observable L0
        error(1.0) D0
        error(1.0) D0 D1 L0
        """))

        @test detector_count(model) == 2
        @test logical_observable_count(model) == 1

        frames = pftrajectories(model; trajectories=5, threads=false)
        @test measurements(frames) == repeat(reshape(Bool[0, 1, 1], 1, 3), 5, 1)
        @test detector_measurements(frames, model) == repeat(reshape(Bool[0, 1], 1, 2), 5, 1)
        @test logical_observable_measurements(frames, model) == trues(5, 1)
    end

    @testset "declarations allocate columns without errors" begin
        declarations_only = detector_error_model_circuit(IOBuffer("""
        detector D3
        logical_observable L1
        """))

        @test detector_count(declarations_only) == 4
        @test logical_observable_count(declarations_only) == 2
        @test isempty(declarations_only)
        @test measurements(pftrajectories(declarations_only; trajectories=3, threads=false)) == falses(3, 6)

        model = detector_error_model_circuit(IOBuffer("""
        detector D3
        logical_observable L1
        error(0.0) D3 L1
        """))

        @test detector_count(model) == 4
        @test logical_observable_count(model) == 2
        @test measurements(pftrajectories(model; trajectories=3, threads=false)) == falses(3, 6)
    end

    @testset "compactification preserves detector model metadata" begin
        model = detector_error_model_circuit(IOBuffer("""
        detector D0
        logical_observable L0
        error(1.0) D0 L0
        """))
        compact_model = compactify_circuit(model)

        @test compact_model isa DetectorErrorModelCircuit
        @test eltype(compact_model) <: QuantumClifford.CompactifiedGate
        @test detector_count(compact_model) == 1
        @test logical_observable_count(compact_model) == 1
        @test measurements(pftrajectories(compact_model; trajectories=3, threads=false)) == trues(3, 2)
    end

    @testset "repeat blocks and detector shifts" begin
        model = detector_error_model_circuit(IOBuffer("""
        repeat 2 {
            repeat 2 {
                error(1.0) D0
                shift_detectors 1
            }
        }
        """))

        @test detector_count(model) == 4
        @test logical_observable_count(model) == 0
        @test measurements(pftrajectories(model; trajectories=4, threads=false)) == trues(4, 4)
    end

    @testset "duplicate targets cancel and separators are ignored" begin
        model = detector_error_model_circuit(IOBuffer("""
        detector D0
        logical_observable L0
        error(1.0) D0 L0 ^ D0 L0
        """))

        @test measurements(pftrajectories(model; trajectories=3, threads=false)) == falses(3, 2)
    end

    @testset "sampled frequencies match error probabilities" begin
        Random.seed!(1234)
        model = detector_error_model_circuit(IOBuffer("error(0.25) D0 L0\n"))
        trajectories = 10_000
        sampled = measurements(pftrajectories(model; trajectories, threads=false))

        @test sampled[:, 1] == sampled[:, 2]
        @test 0.23 < sum(sampled[:, 1]) / trajectories < 0.27
    end

    @testset "path reader and unsupported syntax" begin
        path = tempname()
        write(path, "detector D0\nerror(1.0) D0\n")

        model = read_detector_error_model(path)
        @test detector_count(model) == 1
        @test measurements(pftrajectories(model; trajectories=2, threads=false)) == trues(2, 1)

        @test_throws ArgumentError detector_error_model_circuit(IOBuffer("error(0.1) 5\n"))
    end
end
