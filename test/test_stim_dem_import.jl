@testitem "Stim detector error model import" begin
    using Random

    # Statistical-test design note:
    #
    # Every frequency assertion uses a tolerance of 5 standard errors of the
    # estimated quantity. For a Bernoulli(p) estimate from n trajectories the
    # standard error is SE = sqrt(p*(1-p)/n), so the probability of a false
    # failure is ~6e-7 per assertion for ANY rng seed -- the seeds below are
    # fixed only so that a real failure is reproducible, not to make the tests
    # pass. (For 0/1 data the variance is a deterministic function of the mean,
    # so the informative second-order checks are the cross-column relations,
    # which the models below pin exactly, trajectory by trajectory.)
    se(p, n) = sqrt(p * (1 - p) / n)

    @testset "parsing the minimal example from the issue" begin
        circuit = parse_detector_error_model("""
            detector D0
            detector D1
            logical_observable L0
            error(0.1) D0
            error(0.2) D0 D1 L0
            """)
        @test circuit isa DetectorErrorModelCircuit
        @test circuit isa AbstractVector{QuantumClifford.AbstractOperation}
        @test circuit.n_detectors == 2
        @test circuit.n_logicals == 1
        @test length(circuit) == 3
        @test circuit[1] == DemDeclaration(2, 1)
        # detector Dk -> column k+1; logical Lk -> column n_detectors+k+1
        @test circuit[2].p == 0.1 && circuit[2].detector_bits == [1] && circuit[2].logical_bits == Int[]
        @test circuit[3].p == 0.2 && circuit[3].detector_bits == [1, 2] && circuit[3].logical_bits == [3]
        @test collect(QuantumClifford.affectedbits(circuit[1])) == [1, 2, 3]
        @test QuantumClifford.affectedbits(circuit[3]) == [1, 2, 3]
        @test affectedqubits(circuit[2]) == ()
    end

    @testset "sampling the minimal example: marginals, correlations, output shape" begin
        # Mechanism A ~ Bernoulli(0.1) flips D0; mechanism B ~ Bernoulli(0.2)
        # flips D0, D1, L0. Therefore D0 = A xor B, D1 = L0 = B, and:
        #   P(D0) = 0.1*0.8 + 0.2*0.9 = 0.26
        #   P(D1) = P(L0) = 0.2
        #   P(D0 xor D1) = P(A) = 0.1     (cross-column check, sees correlations)
        #   D1 == L0 exactly, in every single trajectory
        circuit = parse_detector_error_model("""
            detector D0
            detector D1
            logical_observable L0
            error(0.1) D0
            error(0.2) D0 D1 L0
            """)
        n = 20_000
        Random.seed!(728)
        frames = pftrajectories(circuit; trajectories=n, threads=false)
        m = measurements(frames)
        @test size(m) == (n, 3) # trajectories × (n_detectors + n_logicals)
        @test size(detectorview(circuit, frames)) == (n, 2)
        @test size(observableview(circuit, frames)) == (n, 1)
        @test detectorview(circuit, frames) == m[:, 1:2]
        @test observableview(circuit, frames) == m[:, 3:3]
        @test abs(sum(m[:, 1]) / n - 0.26) < 5 * se(0.26, n)  # 5*SE = 0.0155
        @test abs(sum(m[:, 2]) / n - 0.20) < 5 * se(0.20, n)  # 5*SE = 0.0141
        @test abs(sum(m[:, 3]) / n - 0.20) < 5 * se(0.20, n)
        @test abs(sum(m[:, 1] .⊻ m[:, 2]) / n - 0.10) < 5 * se(0.10, n)  # 5*SE = 0.0106
        @test m[:, 2] == m[:, 3] # exact: D1 and L0 are flipped by the same mechanism
        # Same checks through the threaded code path (unseeded; the 5*SE bounds
        # hold for any seed).
        framest = pftrajectories(circuit; trajectories=n, threads=true)
        mt = measurements(framest)
        @test size(mt) == (n, 3)
        @test abs(sum(mt[:, 1]) / n - 0.26) < 5 * se(0.26, n)
        @test mt[:, 2] == mt[:, 3]
    end

    @testset "statistical power guard (falsification check)" begin
        # The bounds above must be able to DETECT a miscalibrated import.
        # Sample a deliberately corrupted version of the model (0.2 -> 0.3);
        # its frequencies must fall OUTSIDE the acceptance bands used above:
        #   P(D0) becomes 0.1*0.7 + 0.3*0.9 = 0.34, which is |0.34-0.26| = 26 SE
        #   away from the correct value, while the band is only 5 SE wide.
        corrupted = parse_detector_error_model("""
            error(0.1) D0
            error(0.3) D0 D1 L0
            """)
        n = 20_000
        Random.seed!(729)
        mc = measurements(pftrajectories(corrupted; trajectories=n, threads=false))
        @test abs(sum(mc[:, 1]) / n - 0.26) > 5 * se(0.26, n)
        @test abs(sum(mc[:, 2]) / n - 0.20) > 5 * se(0.20, n)
    end

    @testset "logical observables land in the columns after the detectors" begin
        # The most likely wiring mistake is an off-by-one (or off-by-n_detectors)
        # in the measurement-column mapping, so pin it deterministically:
        # D = 2 detectors exist, so L0 must land exactly in column 3.
        circuit = parse_detector_error_model("""
            detector D1
            error(1) L0
            """)
        @test circuit.n_detectors == 2 && circuit.n_logicals == 1
        m = measurements(pftrajectories(circuit; trajectories=64, threads=false))
        @test size(m) == (64, 3)
        @test all(.!m[:, 1]) && all(.!m[:, 2]) && all(m[:, 3])
    end

    @testset "repeat blocks and shift_detectors (circular model from the Stim docs)" begin
        # 10 mechanisms with p = 0.1 arranged in a ring over 10 detectors;
        # every detector is flipped by exactly 2 independent mechanisms:
        #   P(detector) = 2 * 0.1 * 0.9 = 0.18, P(L0) = 0.1
        circuit = parse_detector_error_model("""
            error(0.1) D9 D0 L0
            repeat 9 {
                error(0.1) D0 D1
                shift_detectors 1
            }
            """)
        @test circuit.n_detectors == 10
        @test circuit.n_logicals == 1
        @test length(circuit) == 11
        # the unrolled mechanisms walk around the ring
        @test circuit[2].detector_bits == [1, 10]
        @test circuit[3].detector_bits == [1, 2]
        @test circuit[11].detector_bits == [9, 10]
        n = 20_000
        Random.seed!(730)
        m = measurements(pftrajectories(circuit; trajectories=n, threads=false))
        @test size(m) == (n, 11)
        for d in 1:10
            @test abs(sum(m[:, d]) / n - 0.18) < 5 * se(0.18, n)  # 5*SE = 0.0136
        end
        @test abs(sum(m[:, 11]) / n - 0.10) < 5 * se(0.10, n)
    end

    @testset "separators and parity cancellation" begin
        # From the Stim documentation: `error(0.03) D2 L0 ^ D3 L0` has symptoms
        # D2, D3 and NO frame changes (the two L0 cancel out) -- but the
        # cancelled L0 mentions still count towards the model size.
        circuit = parse_detector_error_model("error(0.5) D0 L0 ^ D1 L0")
        @test circuit.n_detectors == 2 && circuit.n_logicals == 1
        @test circuit[2].detector_bits == [1, 2] && circuit[2].logical_bits == Int[]
        m = measurements(pftrajectories(circuit; trajectories=2_000, threads=false))
        @test m[:, 1] == m[:, 2]  # exact: D0 and D1 always flip together
        @test all(.!m[:, 3])      # exact: L0 is never flipped
        @test 0 < sum(m[:, 1]) < 2_000 # and the mechanism does fire sometimes
        # an error whose targets fully cancel flips nothing, but still
        # implicitly declares 4 detectors (Stim's `count_detectors` semantics)
        cancelled = parse_detector_error_model("error(0.5) D3 ^ D3")
        @test cancelled.n_detectors == 4
        @test cancelled[2].detector_bits == Int[]
        mc = measurements(pftrajectories(cancelled; trajectories=128, threads=false))
        @test size(mc) == (128, 4)
        @test !any(mc)
    end

    @testset "declarations allocate output columns; implicit declarations; edge cases" begin
        n = 128
        # declared-but-never-flipped detector must still get a (constant false) column
        c1 = parse_detector_error_model("detector D7\nerror(0.5) D0")
        @test c1.n_detectors == 8
        m1 = measurements(pftrajectories(c1; trajectories=n, threads=false))
        @test size(m1) == (n, 8)
        @test !any(m1[:, 2:8])
        # mentioning D5 in an error implicitly declares D0..D5
        c2 = parse_detector_error_model("error(0.25) D5")
        @test c2.n_detectors == 6 && c2.n_logicals == 0
        @test size(measurements(pftrajectories(c2; trajectories=n, threads=false))) == (n, 6)
        # a trailing shift_detectors does NOT pad the model (matches Stim's
        # count_detectors, where only mentioned/declared indices count)
        c3 = parse_detector_error_model("error(0.5) D0\nshift_detectors 5")
        @test c3.n_detectors == 1
        @test size(measurements(pftrajectories(c3; trajectories=n, threads=false))) == (n, 1)
        # repeat 0 is legal and skips its body (including its shifts)
        c4 = parse_detector_error_model("repeat 0 {\n error(1) D3\n shift_detectors 7\n}\nerror(1) D0")
        @test c4.n_detectors == 1 && length(c4) == 2
        @test all(measurements(pftrajectories(c4; trajectories=n, threads=false))[:, 1])
        # nested repeat blocks with shifts: 2 * 3 iterations, one detector each
        c5 = parse_detector_error_model("""
            repeat 2 {
                repeat 3 {
                    error(1) D0
                    shift_detectors 1
                }
            }
            """)
        @test c5.n_detectors == 6 && length(c5) == 7
        @test [only(c5[i].detector_bits) for i in 2:7] == [1, 2, 3, 4, 5, 6]
        @test all(measurements(pftrajectories(c5; trajectories=n, threads=false)))
        # p = 0 mechanisms never fire
        c6 = parse_detector_error_model("error(0) D0 L0")
        @test !any(measurements(pftrajectories(c6; trajectories=n, threads=false)))
        # empty model
        c7 = parse_detector_error_model("# only a comment\n\n")
        @test c7.n_detectors == 0 && c7.n_logicals == 0 && length(c7) == 1
        @test size(detectorview(c7, pftrajectories(c7; trajectories=n, threads=false))) == (n, 0)
    end

    @testset "tolerated syntax: comments, indentation, tags, case, multi-target declarations" begin
        circuit = parse_detector_error_model("""
            # a comment line
              ERROR[tag with # inside](0.125) D0 d1 L0   # trailing comment
            Detector(1.5, -2, 3) D2 D3
            LOGICAL_OBSERVABLE L1 l2
            Shift_Detectors(0, 0.5) 2
            error(1.0) D1
            REPEAT[t] 2 {
                error(0.25) D0
            }
            """)
        # Mentions: D0, D1 before the shift (absolute 0, 1); declarations D2 D3
        # before the shift (absolute 2, 3); after `shift_detectors 2` the
        # mentions D1 and D0 are absolute 3 and 2. Max absolute index = 3.
        @test circuit.n_detectors == 4
        @test circuit.n_logicals == 3 # L0 from the error, L1 L2 declared
        @test length(circuit) == 5 # declaration + 2 single + 2 unrolled mechanisms
        @test circuit[2].detector_bits == [1, 2] && circuit[2].logical_bits == [5] # L0 -> column 4+0+1
        @test circuit[3].detector_bits == [4] && circuit[3].p == 1.0 # D1 after shift 2 -> absolute 3 -> column 4
        @test circuit[4] == circuit[5] && circuit[4].detector_bits == [3] && circuit[4].p == 0.25 # unrolled repeat
    end

    @testset "clear errors for unsupported or malformed syntax" begin
        @test_throws ArgumentError parse_detector_error_model("frobnicate(0.1) D0")
        @test_throws "unrecognized instruction name 'frobnicate'" parse_detector_error_model("frobnicate(0.1) D0")
        @test_throws "line 2" parse_detector_error_model("error(0.1) D0\nfrobnicate(0.1) D0")
        @test_throws "takes exactly 1 argument" parse_detector_error_model("error D0")
        @test_throws "takes exactly 1 argument" parse_detector_error_model("error(0.1, 0.2) D0")
        @test_throws "must be a probability" parse_detector_error_model("error(1.5) D0")
        @test_throws "must be a probability" parse_detector_error_model("error(-0.5) D0")
        @test_throws "expected a number" parse_detector_error_model("error(zzz) D0")
        @test_throws "unrecognized 'error' target" parse_detector_error_model("error(0.1) D0 X1")
        @test_throws "separators (^)" parse_detector_error_model("error(0.1) ^ D0")
        @test_throws "separators (^)" parse_detector_error_model("error(0.1) D0 ^")
        @test_throws "separators (^)" parse_detector_error_model("error(0.1) D0 ^ ^ D1")
        @test_throws "expected a non-negative integer" parse_detector_error_model("error(0.1) D-1")
        @test_throws "relative detector targets" parse_detector_error_model("detector L0")
        @test_throws "at least 1 target" parse_detector_error_model("detector")
        @test_throws "takes 0 arguments" parse_detector_error_model("logical_observable(1) L0")
        @test_throws "logical observable targets" parse_detector_error_model("logical_observable D0")
        @test_throws "exactly 1 numeric target" parse_detector_error_model("shift_detectors")
        @test_throws "exactly 1 numeric target" parse_detector_error_model("shift_detectors 1 2")
        @test_throws "expected a non-negative integer" parse_detector_error_model("shift_detectors D1")
        @test_throws "missing '{'" parse_detector_error_model("repeat 5")
        @test_throws "unterminated block" parse_detector_error_model("repeat 5 {\nerror(0.1) D0")
        @test_throws "got a '}' without a '{'" parse_detector_error_model("}")
        @test_throws "unexpected '{'" parse_detector_error_model("error(0.1) D0 {\n}")
        @test_throws "unterminated instruction tag" parse_detector_error_model("error[oops(0.1) D0")
        # constructor validation
        @test_throws ArgumentError DetectorError(1.5, [1], Int[])
        @test_throws ArgumentError DetectorError(0.5, [0], Int[])
        @test_throws ArgumentError DemDeclaration(-1, 0)
    end

    @testset "reading from a file and from an IO stream" begin
        content = """
            error(0.1) D0
            error(0.2) D0 D1 L0
            """
        mktempdir() do dir
            path = joinpath(dir, "model.dem")
            write(path, content)
            circuit = read_detector_error_model(path)
            @test circuit isa DetectorErrorModelCircuit
            @test circuit.n_detectors == 2 && circuit.n_logicals == 1
            @test circuit == parse_detector_error_model(content)
            open(path) do io
                @test read_detector_error_model(io) == circuit
            end
        end
    end

    @testset "cross-check against mctrajectories on a Register" begin
        # The same DetectorError operations also act on the classical bits of a
        # Register, so the Monte Carlo backend provides an independent check of
        # the Pauli-frame sampling statistics (and of `affectedbits` wiring).
        circuit = parse_detector_error_model("""
            error(0.1) D0
            error(0.2) D0 D1 L0
            """)
        n = 4_000
        Random.seed!(731)
        reg = Register(one(MixedDestabilizer, 1), 3)
        bits = stack([bitview(mctrajectory!(copy(reg), circuit)[1]) for _ in 1:n], dims=1)
        @test size(bits) == (n, 3)
        @test abs(sum(bits[:, 1]) / n - 0.26) < 5 * se(0.26, n)
        @test abs(sum(bits[:, 2]) / n - 0.20) < 5 * se(0.20, n)
        @test bits[:, 2] == bits[:, 3]
    end

    @testset "a realistic folded-loop file (Stim docs repetition code, d=4, r=1000)" begin
        # Abridged structural copy of the `--gen repetition_code --rounds 1000
        # --distance 4` example from Stim's .dem format documentation
        # (probabilities rounded; coordinates kept). It exercises coordinate
        # arguments, interleaved shifts inside a 998-iteration loop, and
        # trailing declarations. Expected: 3 detectors per round over 1001
        # rounds = 3003 detectors, 13 mechanisms per round block = 13000, plus
        # the leading declaration operation.
        circuit = parse_detector_error_model("""
            error(0.00027) D0
            error(0.00027) D0 D1
            error(0.00053) D0 D3
            error(0.00053) D0 D4
            error(0.00027) D1 D2
            error(0.00053) D1 D4
            error(0.00053) D1 D5
            error(0.00053) D2 D5
            error(0.00027) D2 L0
            error(0.00027) D3
            error(0.00027) D3 D4
            error(0.00027) D4 D5
            error(0.00027) D5 L0
            detector(1, 0) D0
            detector(3, 0) D1
            detector(5, 0) D2
            repeat 998 {
                error(0.00027) D3
                error(0.00027) D3 D4
                error(0.00053) D3 D6
                error(0.00053) D3 D7
                error(0.00027) D4 D5
                error(0.00053) D4 D7
                error(0.00053) D4 D8
                error(0.00053) D5 D8
                error(0.00027) D5 L0
                error(0.00027) D6
                error(0.00027) D6 D7
                error(0.00027) D7 D8
                error(0.00027) D8 L0
                shift_detectors(0, 1) 0
                detector(1, 0) D3
                detector(3, 0) D4
                detector(5, 0) D5
                shift_detectors 3
            }
            error(0.00027) D3
            error(0.00027) D3 D4
            error(0.00053) D3 D6
            error(0.00053) D3 D7
            error(0.00027) D4 D5
            error(0.00053) D4 D7
            error(0.00053) D4 D8
            error(0.00053) D5 D8
            error(0.00027) D5 L0
            error(0.00027) D6
            error(0.00027) D6 D7
            error(0.00027) D7 D8
            error(0.00027) D8 L0
            shift_detectors(0, 1) 0
            detector(1, 0) D3
            detector(3, 0) D4
            detector(5, 0) D5
            detector(1, 1) D6
            detector(3, 1) D7
            detector(5, 1) D8
            """)
        @test circuit.n_detectors == 3003
        @test circuit.n_logicals == 1
        @test length(circuit) == 13_001
        # the last unrolled mechanism is `error D8 L0` of the final round,
        # at detector offset 2994: absolute D3002 -> column 3003
        @test circuit[end].detector_bits == [3003]
        @test circuit[end].logical_bits == [3004]
        frames = pftrajectories(circuit; trajectories=50, threads=false)
        @test size(measurements(frames)) == (50, 3004)
    end
end
