@testitem "Detector error model (.dem) import" begin
    using QuantumClifford
    using QuantumClifford: DetectorError, DetectorErrorModelColumns, DetectorErrorModelCircuit

    read_dem(str) = read_detector_error_model(IOBuffer(str))

    @testset "minimal example from the issue" begin
        dem = """
        detector D0
        detector D1
        logical_observable L0
        error(0.1) D0
        error(0.2) D0 D1 L0
        """
        circuit = read_dem(dem)
        @test circuit isa DetectorErrorModelCircuit
        @test circuit.num_detectors == 2
        @test circuit.num_logicals == 1

        n = 2_000_000
        frames = pftrajectories(circuit; trajectories=n)
        m = measurements(frames)
        @test size(m) == (n, 3) # 2 detectors + 1 logical

        D0 = m[:, 1]; D1 = m[:, 2]; L0 = m[:, 3]
        # D0 is flipped by both mechanisms (p=0.1 and p=0.2, independently)
        pD0 = 0.1*(1-0.2) + 0.2*(1-0.1)
        # D1 and L0 are flipped only by the second mechanism (p=0.2)
        @test isapprox(sum(D0)/n, pD0; atol=0.005)
        @test isapprox(sum(D1)/n, 0.2; atol=0.005)
        @test isapprox(sum(L0)/n, 0.2; atol=0.005)
        # D1 and L0 are perfectly correlated (always flipped together)
        @test all(D1 .== L0)
    end

    @testset "single independent mechanism statistics" begin
        for p in (0.05, 0.3, 0.5, 0.9)
            circuit = read_dem("error($p) D0")
            n = 1_000_000
            m = measurements(pftrajectories(circuit; trajectories=n))
            @test size(m, 2) == 1
            @test isapprox(sum(m[:,1])/n, p; atol=0.005)
        end
    end

    @testset "two mechanisms on the same detector combine correctly" begin
        # p_combined = p1(1-p2) + p2(1-p1)
        p1, p2 = 0.1, 0.25
        circuit = read_dem("error($p1) D0\nerror($p2) D0")
        n = 2_000_000
        m = measurements(pftrajectories(circuit; trajectories=n))
        expected = p1*(1-p2) + p2*(1-p1)
        @test isapprox(sum(m[:,1])/n, expected; atol=0.005)
    end

    @testset "separators only suggest a decomposition (parity semantics)" begin
        # error(p) D2 L0 ^ D3 L0 -> symptoms D2, D3, no frame change (the two L0 cancel)
        circuit = read_dem("error(0.3) D2 L0 ^ D3 L0")
        @test circuit.num_detectors == 4   # D0..D3
        @test circuit.num_logicals == 1    # L0 mentioned (even if it cancels)
        n = 1_000_000
        m = measurements(pftrajectories(circuit; trajectories=n))
        D2 = m[:, 3]; D3 = m[:, 4]; L0 = m[:, 5]
        @test isapprox(sum(D2)/n, 0.3; atol=0.005)
        @test all(D2 .== D3)          # always flipped together
        @test sum(L0) == 0            # L0 never flips (cancelled out)
    end

    @testset "repeat blocks with shift_detectors (circular model)" begin
        # Equivalent to the circular error model from the Stim docs.
        dem_repeat = """
        error(0.1) D9 D0 L0
        repeat 9 {
            error(0.1) D0 D1
            shift_detectors 1
        }
        """
        dem_flat = """
        error(0.1) D9 D0 L0
        error(0.1) D0 D1
        error(0.1) D1 D2
        error(0.1) D2 D3
        error(0.1) D3 D4
        error(0.1) D4 D5
        error(0.1) D5 D6
        error(0.1) D6 D7
        error(0.1) D7 D8
        error(0.1) D8 D9
        """
        cr = read_dem(dem_repeat)
        cf = read_dem(dem_flat)
        @test cr.num_detectors == 10 == cf.num_detectors
        @test cr.num_logicals == 1 == cf.num_logicals
        # both should produce 10 error mechanisms
        @test count(o->isa(o, DetectorError), cr.ops) == 10
        @test count(o->isa(o, DetectorError), cf.ops) == 10
    end

    @testset "nested repeat blocks" begin
        dem = """
        repeat 2 {
            repeat 3 {
                error(0.1) D0
                shift_detectors 1
            }
        }
        """
        circuit = read_dem(dem)
        # 2*3 = 6 mechanisms, detectors D0..D5
        @test count(o->isa(o, DetectorError), circuit.ops) == 6
        @test circuit.num_detectors == 6
    end

    @testset "comments, blank lines, tags, and coordinates are tolerated" begin
        dem = """
        # a leading comment
        detector(1, 0) D0      # detector with coordinates

        error[some_tag](0.2) D0   # an error with a tag and a comment
        shift_detectors(0, 0.5) 0 # only coordinate shift, no detector shift
        logical_observable L0
        """
        circuit = read_dem(dem)
        @test circuit.num_detectors == 1
        @test circuit.num_logicals == 1
        @test count(o->isa(o, DetectorError), circuit.ops) == 1
    end

    @testset "declared-but-unused detectors/observables still allocate columns" begin
        dem = """
        detector D4
        logical_observable L2
        error(0.5) D0
        """
        circuit = read_dem(dem)
        @test circuit.num_detectors == 5  # D0..D4
        @test circuit.num_logicals == 3   # L0..L2
        m = measurements(pftrajectories(circuit; trajectories=1000))
        @test size(m, 2) == 8             # 5 detectors + 3 logicals
        # only D0 ever flips
        @test all(sum(m, dims=1)[1, 2:end] .== 0)
    end

    @testset "case-insensitive instruction names and targets" begin
        circuit = read_dem("ERROR(0.5) d0 l0")
        @test circuit.num_detectors == 1
        @test circuit.num_logicals == 1
        @test count(o->isa(o, DetectorError), circuit.ops) == 1
    end

    @testset "clear errors for unsupported / malformed syntax" begin
        @test_throws ArgumentError read_dem("not_a_real_instruction D0")
        @test_throws ArgumentError read_dem("error(2.0) D0")        # probability out of range
        @test_throws ArgumentError read_dem("error D0")             # missing probability
        @test_throws ArgumentError read_dem("error(0.1) Q0")        # invalid target
        @test_throws ArgumentError read_dem("repeat 3 {\nerror(0.1) D0")  # missing closing brace
        @test_throws ArgumentError read_dem("}")                    # stray closing brace
    end

    @testset "reading from a file path" begin
        dem = "error(0.5) D0 L0\n"
        path, io = mktemp()
        write(io, dem); close(io)
        circuit = read_detector_error_model(path)
        @test circuit.num_detectors == 1
        @test circuit.num_logicals == 1
        rm(path; force=true)
    end
end
