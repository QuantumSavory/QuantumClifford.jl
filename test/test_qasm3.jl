@testitem "OpenQASM 3 import (QuantumCliffordQuasarExt)" tags=[:quasar] begin
    # This test file requires Quasar.jl to be available.
    # It is gated so that the CI matrix can skip it when Quasar is not installed.
    using QuantumClifford
    import Quasar  # triggers QuantumCliffordQuasarExt

    FIXTURES = joinpath(@__DIR__, "qasm_fixtures")

    # -----------------------------------------------------------------------
    # 1.  Bell circuit: the canonical example from the issue
    # -----------------------------------------------------------------------
    @testset "Bell circuit from string" begin
        src = """
        OPENQASM 3.0;
        qubit[2] q;
        bit[2] c;
        h q[0];
        cx q[0], q[1];
        c[0] = measure q[0];
        c[1] = measure q[1];
        """
        circuit = parse_qasm3(src)

        @test length(circuit) == 4
        @test circuit[1] == sHadamard(1)
        @test circuit[2] == sCNOT(1, 2)
        @test circuit[3] == sMZ(1, 1)
        @test circuit[4] == sMZ(2, 2)
    end

    @testset "Bell circuit from file" begin
        circuit = read_qasm3(joinpath(FIXTURES, "bell.qasm"))
        @test length(circuit) == 4
        @test circuit[1] == sHadamard(1)
        @test circuit[2] == sCNOT(1, 2)
        @test circuit[3] == sMZ(1, 1)
        @test circuit[4] == sMZ(2, 2)
    end

    # -----------------------------------------------------------------------
    # 2.  Bell circuit runs through pftrajectories without error
    # -----------------------------------------------------------------------
    @testset "Bell circuit pftrajectories" begin
        circuit = read_qasm3(joinpath(FIXTURES, "bell.qasm"))
        @test_nowarn pftrajectories(circuit; trajectories=100)
    end

    # -----------------------------------------------------------------------
    # 3.  Single-qubit gate mapping
    # -----------------------------------------------------------------------
    @testset "Single-qubit gate mapping" begin
        single_q_map = [
            ("id q[0];", sId1(1)),
            ("x q[0];",  sX(1)),
            ("y q[0];",  sY(1)),
            ("z q[0];",  sZ(1)),
            ("h q[0];",  sHadamard(1)),
            ("s q[0];",  sPhase(1)),
            ("sdg q[0];", sInvPhase(1)),
        ]
        for (gate_line, expected) in single_q_map
            src = """
            OPENQASM 3.0;
            qubit[1] q;
            $gate_line
            """
            circuit = parse_qasm3(src)
            @test length(circuit) == 1
            @test circuit[1] == expected
        end
    end

    # -----------------------------------------------------------------------
    # 4.  Two-qubit gate mapping
    # -----------------------------------------------------------------------
    @testset "Two-qubit gate mapping" begin
        two_q_map = [
            ("cx   q[0], q[1];",  sCNOT(1, 2)),
            ("cz   q[0], q[1];",  sCPHASE(1, 2)),
            ("swap q[0], q[1];",  sSWAP(1, 2)),
        ]
        for (gate_line, expected) in two_q_map
            src = """
            OPENQASM 3.0;
            qubit[2] q;
            $gate_line
            """
            circuit = parse_qasm3(src)
            @test length(circuit) == 1
            @test circuit[1] == expected
        end
    end

    # -----------------------------------------------------------------------
    # 5.  Classical bit array and qubit array declarations
    # -----------------------------------------------------------------------
    @testset "Qubit and bit array declarations" begin
        src = """
        OPENQASM 3.0;
        qubit[3] q;
        bit[3]   c;
        h q[0];
        x q[1];
        z q[2];
        c[0] = measure q[0];
        c[1] = measure q[1];
        c[2] = measure q[2];
        """
        circuit = parse_qasm3(src)
        @test length(circuit) == 6
        @test circuit[1] == sHadamard(1)
        @test circuit[2] == sX(2)
        @test circuit[3] == sZ(3)
        @test circuit[4] == sMZ(1, 1)
        @test circuit[5] == sMZ(2, 2)
        @test circuit[6] == sMZ(3, 3)
    end

    # -----------------------------------------------------------------------
    # 6.  Measurement assignment: correct qubit → bit index binding
    # -----------------------------------------------------------------------
    @testset "Measurement assignment" begin
        src = """
        OPENQASM 3.0;
        qubit[2] q;
        bit[2]   c;
        c[1] = measure q[0];
        c[0] = measure q[1];
        """
        circuit = parse_qasm3(src)
        @test length(circuit) == 2
        @test circuit[1] == sMZ(1, 2)   # qubit 1 → classical bit 2
        @test circuit[2] == sMZ(2, 1)   # qubit 2 → classical bit 1
    end

    # -----------------------------------------------------------------------
    # 7.  Reset operation
    # -----------------------------------------------------------------------
    @testset "Reset operation" begin
        src = """
        OPENQASM 3.0;
        qubit[1] q;
        reset q[0];
        """
        circuit = parse_qasm3(src)
        @test length(circuit) == 1
        @test circuit[1] == sMRZ(1, 0)
    end

    @testset "Reset circuit from file" begin
        circuit = read_qasm3(joinpath(FIXTURES, "reset.qasm"))
        # h, cx, reset q[0], reset q[1], measure q[0], measure q[1]
        @test length(circuit) == 6
        @test circuit[3] == sMRZ(1, 0)
        @test circuit[4] == sMRZ(2, 0)
    end

    # -----------------------------------------------------------------------
    # 8.  Multiple registers — correct offset arithmetic
    # -----------------------------------------------------------------------
    @testset "Multiple qubit and bit registers" begin
        circuit = read_qasm3(joinpath(FIXTURES, "multi_register.qasm"))
        # qubit[2] a  → qubits 1,2
        # qubit[2] b  → qubits 3,4
        # bit[2]  ca  → bits  1,2
        # bit[2]  cb  → bits  3,4
        # h a[0]        → sHadamard(1)
        # cx a[0], b[0] → sCNOT(1,3)
        # x b[1]        → sX(4)
        # ca[0] = measure a[0] → sMZ(1,1)
        # ca[1] = measure a[1] → sMZ(2,2)
        # cb[0] = measure b[0] → sMZ(3,3)
        # cb[1] = measure b[1] → sMZ(4,4)
        @test length(circuit) == 7
        @test circuit[1] == sHadamard(1)
        @test circuit[2] == sCNOT(1, 3)
        @test circuit[3] == sX(4)
        @test circuit[4] == sMZ(1, 1)
        @test circuit[5] == sMZ(2, 2)
        @test circuit[6] == sMZ(3, 3)
        @test circuit[7] == sMZ(4, 4)
    end

    # -----------------------------------------------------------------------
    # 9.  0-based → 1-based index conversion
    # -----------------------------------------------------------------------
    @testset "Zero-based to one-based conversion" begin
        # q[0] in QASM should become qubit 1 in QuantumClifford
        src = """
        OPENQASM 3.0;
        qubit[4] q;
        x q[0];
        x q[3];
        """
        circuit = parse_qasm3(src)
        @test circuit[1] == sX(1)
        @test circuit[2] == sX(4)
    end

    # -----------------------------------------------------------------------
    # 10.  All supported gates file
    # -----------------------------------------------------------------------
    @testset "All Clifford gates fixture" begin
        circuit = read_qasm3(joinpath(FIXTURES, "all_clifford_gates.qasm"))
        # Should parse without errors; verify length
        # id, x, y, z, h, s, sdg (7) + cx, cz, swap (3) + 6 measurements = 16
        @test length(circuit) == 16
    end

    # -----------------------------------------------------------------------
    # 11.  Unsupported gate → informative error
    # -----------------------------------------------------------------------
    @testset "Unsupported gate error" begin
        src = """
        OPENQASM 3.0;
        qubit[1] q;
        t q[0];
        """
        err = @test_throws ArgumentError parse_qasm3(src)
        @test occursin("t", lowercase(err.value.msg))
    end

    @testset "Unsupported gate from file" begin
        err = @test_throws ArgumentError read_qasm3(
            joinpath(FIXTURES, "unsupported_t_gate.qasm")
        )
        @test occursin("t", lowercase(err.value.msg))
    end

    @testset "Unsupported parameterized gate error" begin
        src = """
        OPENQASM 3.0;
        qubit[1] q;
        rx(1.5708) q[0];
        """
        @test_throws ArgumentError parse_qasm3(src)
    end

    # -----------------------------------------------------------------------
    # 12.  Round-trip: imported Bell circuit simulation produces correlations
    # -----------------------------------------------------------------------
    @testset "Bell simulation correlations via pftrajectories" begin
        circuit = parse_qasm3("""
        OPENQASM 3.0;
        qubit[2] q;
        bit[2]   c;
        h  q[0];
        cx q[0], q[1];
        c[0] = measure q[0];
        c[1] = measure q[1];
        """)
        N = 10_000
        frames = pftrajectories(circuit; trajectories=N)
        meas = pfmeasurements(frames)    # N×2 BitMatrix
        # In a Bell state: both bits should be equal (00 or 11), never 01 or 10
        @test all(meas[:, 1] .== meas[:, 2])
    end

end
