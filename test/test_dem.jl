@testitem "Detector error model import" begin
    using QuantumClifford
    using QuantumClifford: DetectorError, DeclareMeasurementBits
    using Random

    @testset "parsing the minimal example" begin
        dem = """
        detector D0
        detector D1
        logical_observable L0
        error(0.1) D0
        error(0.2) D0 D1 L0
        """
        c = parse_detector_error_model(dem)
        @test c isa DetectorErrorModelCircuit
        @test c.num_detectors == 2
        @test c.num_observables == 1
        errs = [op for op in c if op isa DetectorError]
        @test length(errs) == 2
        # detectors map to columns 1:2, the single observable to column 3
        @test errs[1].p == 0.1 && errs[1].bits == (1,)
        @test errs[2].p == 0.2 && errs[2].bits == (1, 2, 3)
        # a declaration is appended so the buffer is always the full width
        @test any(op isa DeclareMeasurementBits && op.nbits == 3 for op in c)
    end

    @testset "comments, blank lines, tags and decomposition separators" begin
        dem = """
        # leading comment

        error[a_tag](0.3) D0 ^ D1 L0  # trailing comment
        logical_observable L2
        """
        c = parse_detector_error_model(dem)
        @test c.num_detectors == 2   # D0, D1
        @test c.num_observables == 3 # L0, L1, L2 -> highest index 2
        err = only(op for op in c if op isa DetectorError)
        # the `^` separator is ignored; the mechanism flips all listed targets
        @test err.bits == (1, 2, 3) # D0->1, D1->2, L0-> num_detectors+0+1 = 3
    end

    @testset "shift_detectors and repeat expansion" begin
        dem = """
        repeat 3 {
            error(0.05) D0 D1
            shift_detectors 1
        }
        """
        c = parse_detector_error_model(dem)
        @test c.num_detectors == 4 # detectors 0..3 across the three shifted rounds
        cols = [op.bits for op in c if op isa DetectorError]
        @test cols == [(1, 2), (2, 3), (3, 4)]
    end

    @testset "nested repeat" begin
        dem = """
        repeat 2 {
          repeat 2 {
            error(0.1) D0
            shift_detectors 1
          }
        }
        """
        c = parse_detector_error_model(dem)
        cols = [op.bits for op in c if op isa DetectorError]
        @test cols == [(1,), (2,), (3,), (4,)]
        @test c.num_detectors == 4
    end

    @testset "declared-but-unused detectors still allocate columns" begin
        dem = """
        detector D3
        error(0.5) D0
        """
        c = parse_detector_error_model(dem)
        @test c.num_detectors == 4 # D0..D3, even though only D0 ever flips
        frames = pftrajectories(c; trajectories=1000, threads=false)
        m = measurements(frames)
        @test size(m) == (1000, 4)
        @test all(.!m[:, 2]) && all(.!m[:, 3]) && all(.!m[:, 4]) # never flipped
    end

    @testset "file round-trip" begin
        path = tempname() * ".dem"
        write(path, "detector D0\nerror(0.5) D0\n")
        @test read_detector_error_model(path) isa DetectorErrorModelCircuit
        @test detector_error_model_circuit(path) isa DetectorErrorModelCircuit
    end

    @testset "clear errors for unsupported syntax" begin
        @test_throws ArgumentError parse_detector_error_model("error(0.1) D0\n}")          # unmatched }
        @test_throws ArgumentError parse_detector_error_model("repeat 2 {\nerror(0.1) D0") # unmatched {
        @test_throws ArgumentError parse_detector_error_model("frobnicate D0")             # unknown instruction
        @test_throws ArgumentError parse_detector_error_model("error(0.1) X3")             # bad target
        @test_throws ArgumentError parse_detector_error_model("error D0")                  # missing probability
        @test_throws ArgumentError parse_detector_error_model("error(1.5) D0")             # probability out of range
        @test_throws ArgumentError parse_detector_error_model("# only comments")           # nothing declared
    end

    @testset "statistical validation against analytic distributions" begin
        Random.seed!(1234)
        n = 2 * 10^5

        # Two independent mechanisms over D0, D1 and L0.
        dem = """
        detector D0
        detector D1
        logical_observable L0
        error(0.1) D0
        error(0.2) D0 D1 L0
        """
        c = parse_detector_error_model(dem)
        for threads in (false, true)
            m = measurements(pftrajectories(c; trajectories=n, threads))
            @test size(m) == (n, 3)
            freqs = vec(sum(m, dims=1)) ./ n
            # D0 flips if exactly one of the two mechanisms fires: 0.1*0.8 + 0.9*0.2
            # D1 and L0 flip only with the second mechanism: 0.2
            expected = [0.1 * 0.8 + 0.9 * 0.2, 0.2, 0.2]
            @test all(abs.(freqs .- expected) .< 0.01)
            # D1 and L0 always flip together (same mechanism)
            @test m[:, 2] == m[:, 3]
        end

        # Correlated detectors: a single mechanism flips D0 and D1 jointly.
        c2 = parse_detector_error_model("error(0.3) D0 D1\n")
        m2 = measurements(pftrajectories(c2; trajectories=n, threads=false))
        @test m2[:, 1] == m2[:, 2]                       # perfectly correlated
        @test abs(sum(m2[:, 1]) / n - 0.3) < 0.01        # correct marginal
    end
end
