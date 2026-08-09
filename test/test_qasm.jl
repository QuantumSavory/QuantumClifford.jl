@testitem "Read OpenQASM" begin
    using Quasar: QasmVisitorError
    using QuantumClifford

    @testset "bell circuit" begin
        @test read_qasm3("test_qasm/test_qasm1.qasm") == [
            sHadamard(1),   # h q[0]; 
            sCNOT(1, 2),    # cx q[0], q[1];
            sMZ(1, 1),      # c[0] = measure q[0]; 
            sMZ(2, 2)       # c[1] = measure q[1];
        ]
    end

    @testset "all supported operations" begin
        @test read_qasm3("test_qasm/test_qasm2.qasm") == [
            # Single qubit tests
            sHadamard(1),                                           # h q1[0]
            sHadamard(5), sHadamard(6), sHadamard(7), sHadamard(8), # h q2;
            sHadamard(4), sHadamard(3),                             # h q1[{3,2}]
            sHadamard(1), sHadamard(2), sHadamard(3),               # h q1[0:2]
            sHadamard(4),                                           # h q1[-1]

            # Two qubit tests
            sCNOT(1,2),                                         # cx q1[0], q1[1]
            sCNOT(1,5), sCNOT(2,6), sCNOT(3,7), sCNOT(4,8),     # cx q1, q2;
            sCNOT(5,1), sCNOT(5,2), sCNOT(5,3), sCNOT(5,4),     # cx q2[0], q1;
            sCNOT(1,5), sCNOT(2,5), sCNOT(3,5), sCNOT(4,5),     # cx q1, q2[0];
            sCNOT(1,2), sCNOT(3,4),                             # cx q1[{0,2}], q1[{1,3}]
            sCNOT(1,2), sCNOT(2,3), sCNOT(3,4),                 # cx q1[0:2], q1[1:3]
            sCNOT(4, 7),                                        # cx q1[-1], q2[-2]
            
            # Test all gates
            sId1(1),        # id q1[0]
            sX(2),          # x q1[1]
            sY(3),          # y q1[2]
            sZ(4),          # z q1[3]
            sPhase(5),      # s q2[0]
            sInvPhase(6),   # sdg q2[5]

            sCPHASE(3,2),   # cz q1[2], q1[1]
            sSWAP(7,6),     # swap q2[2], q2[1]
            
            # Test measurements and resets
            sMZ(4),                                         # measure q1[3]                  
            sMZ(1, 5),                                      # c2 = measure q1[0]
            sMZ(1, 1), sMZ(2, 2), sMZ(3, 3), sMZ(4, 4),     # c1 = measure q1
            sMZ(2, 2), sMZ(3, 3), sMZ(4, 4),                # c1[1:3] = measure q1[1:3]
            sMZ(8, 1), sMZ(7, 3),                           # c1[{0,2}] = measure q2[{3,2}]
            sMZ(2, 3),                                      # c1[-2] = measure q1[-3];
            sMZ(7, 6),                                      # bit c = measure q2[2];
            
            sMRZ(1, 0),                                     # reset q1[0]
            sMRZ(1, 0), sMRZ(2, 0), sMRZ(3, 0), sMRZ(4, 0), # reset q1
            sMRZ(2, 0), sMRZ(3, 0), sMRZ(4, 0),             # reset q1[1:3]
            sMRZ(8, 0), sMRZ(7, 0)                          # reset q2[{3,2}]
        ]                         
    end

    @testset "invalid operations" begin
        # Indexing errors
        @test_throws QasmVisitorError parse_qasm3("""
            bit[2] c;
            qubit[2] q;
            c[2] = measure q[0];
        """)
        @test_throws QasmVisitorError parse_qasm3("""
            bit[2] c;
            qubit[2] q;
            c[-5] = measure q[0];
        """)
        @test_throws QasmVisitorError parse_qasm3("""
            bit[2] c;
            qubit[2] q;
            c[0] = measure q[4];
        """)
        @test_throws QasmVisitorError parse_qasm3("""
            bit[2] c;
            qubit[2] q;
            c[1] = measure q[-5];
        """)

        # Measurement source/destination size mismatch
        @test_throws QasmVisitorError parse_qasm3("""
            qubit[2] q;
            bit[1] c;
            c = measure q;
        """)
        @test_throws QasmVisitorError parse_qasm3("""
            bit[2] c;
            qubit[2] q;
            c[0] = measure q[0:5];
        """)
    end

    @testset "unsupported operations" begin
        # Wrong OpenQASM version
        @test_throws QasmVisitorError parse_qasm3("OPENQASM 2.0;")

        # Unsupported gate
        @test_throws QasmVisitorError parse_qasm3("""
            qubit[1] q;
            t q[0];
        """)

        # Unsupported parameterized gate
        @test_throws QasmVisitorError parse_qasm3("""
            qubit[1] q;
            rx(pi/2) q[0];
        """)

        # Unsupported classical datatype
        @test_throws QasmVisitorError parse_qasm3("bool x;")
        @test_throws QasmVisitorError parse_qasm3("int x;")
        @test_throws QasmVisitorError parse_qasm3("int[16] x;")
        @test_throws QasmVisitorError parse_qasm3("float[64] x;")
        @test_throws QasmVisitorError parse_qasm3("angle[32] x;")
        @test_throws QasmVisitorError parse_qasm3("complex[float[64]] x;")
        @test_throws QasmVisitorError parse_qasm3("duration x;")
        @test_throws QasmVisitorError parse_qasm3("stretch x;")

        # Unsupported classical assignment
        @test_throws QasmVisitorError parse_qasm3("""
            bit[1] c;
            c[0] = true;
        """)
        @test_throws QasmVisitorError parse_qasm3("""bit[2] c = "00";""")

        # Unsupported classical expressions, arithmetic, and control flow
        @test_throws QasmVisitorError parse_qasm3("""
            bit[1] c;
            c[0] = !c[0];
        """)
        @test_throws QasmVisitorError parse_qasm3("""
            bit[2] c;
            c[0] = c[0] ^ c[1];
        """)
        @test_throws QasmVisitorError parse_qasm3("""
            bit[2] c;
            c[0] = c[0] + c[1];
        """)
        @test_throws QasmVisitorError parse_qasm3("""
            qubit[1] q;
            bit[1] c;
            if (c[0]) {
                h q[0];
            }
        """)
        @test_throws QasmVisitorError parse_qasm3("""
            qubit[1] q;
            bit[1] c;
            while (c[0]) {
                c[0] = measure q[0];
            }
        """)

        # Unsupported timing, calibration, and hardware-specific constructs
        @test_throws QasmVisitorError parse_qasm3("""
            qubit[1] q;
            delay[10ns] q;
        """)
        @test_throws QasmVisitorError parse_qasm3("defcal x \$0 {}")

        # Unsupported external file dependencies
        @test_throws SystemError parse_qasm3("""include "stdgates.inc";""")
    end
end