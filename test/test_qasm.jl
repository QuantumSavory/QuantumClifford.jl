@testitem "OpenQASM 3 import" begin
    using QuantumClifford
    import Quasar

    bell_qasm = """
    OPENQASM 3.0;
    qubit[2] q;
    bit[2] c;

    h q[0];
    cx q[0], q[1];
    c[0] = measure q[0];
    c[1] = measure q[1];
    """

    bell_circuit = parse_qasm3(bell_qasm)
    @test bell_circuit == QuantumClifford.AbstractOperation[
        sHadamard(1),
        sCNOT(1, 2),
        sMZ(1, 1),
        sMZ(2, 2),
    ]

    frames = pftrajectories(bell_circuit; trajectories=8, threads=false)
    @test size(measurements(frames)) == (8, 2)

    mktemp() do path, io
        write(io, bell_qasm)
        close(io)
        @test read_qasm3(path) == bell_circuit
    end

    all_supported_qasm = """
    OPENQASM 3.0;
    qubit[3] q;
    bit[3] c;

    id q[0];
    x q[0];
    y q[1];
    z q[2];
    h q[0];
    s q[1];
    sdg q[2];
    cx q[0], q[1];
    cz q[1], q[2];
    swap q[0], q[2];
    reset q[2];
    c[0:2] = measure q[0:2];
    """

    all_supported_circuit = @test_logs (:warn, r"reset expression encountered") parse_qasm3(all_supported_qasm)
    @test all_supported_circuit == QuantumClifford.AbstractOperation[
        sId1(1),
        sX(1),
        sY(2),
        sZ(3),
        sHadamard(1),
        sPhase(2),
        sInvPhase(3),
        sCNOT(1, 2),
        sCPHASE(2, 3),
        sSWAP(1, 3),
        sMRZ(3),
        sMZ(1, 1),
        sMZ(2, 2),
        sMZ(3, 3),
    ]

    custom_gate_qasm = """
    OPENQASM 3.0;
    gate bell a, b {
        h a;
        cx a, b;
    }
    qubit[2] q;
    bell q[0], q[1];
    measure q[0];
    """

    @test parse_qasm3(custom_gate_qasm) == QuantumClifford.AbstractOperation[
        sHadamard(1),
        sCNOT(1, 2),
        sMZ(1),
    ]

    unsupported_gate_err = try
        parse_qasm3("""
        OPENQASM 3.0;
        qubit q;
        t q;
        """)
    catch err
        err
    end
    @test unsupported_gate_err isa Quasar.QasmVisitorError
    @test occursin("unsupported or undefined OpenQASM gate `t`", sprint(showerror, unsupported_gate_err))
    @test occursin("id, x, y, z, h, s, sdg, cx, cz, swap, measure, reset", sprint(showerror, unsupported_gate_err))

    mismatch_err = try
        parse_qasm3("""
        OPENQASM 3.0;
        qubit[2] q;
        bit[1] c;
        c[0] = measure q[0:1];
        """)
    catch err
        err
    end
    @test mismatch_err isa Quasar.QasmVisitorError
    @test occursin("measurement source and destination sizes must match", sprint(showerror, mismatch_err))
end
