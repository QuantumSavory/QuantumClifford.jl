@testitem "Symbolic Clifford" begin
    using Random
    using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
    using QuantumClifford: apply_single_x!, apply_single_y!, apply_single_z!
    using InteractiveUtils
    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    @testset "Small symbolic operators" begin
        for n in test_sizes
            for i in 1:6
                op = enumerate_single_qubit_gates(i, qubit=n, phases=(rand(Bool),rand(Bool)))
                op0 = enumerate_single_qubit_gates(i, qubit=n)
                op_cc = CliffordOperator(op, 1, compact=true)
                op_c = CliffordOperator(op, n)
                @test SingleQubitOperator(op)==SingleQubitOperator(op_cc, n)
                op0_c = CliffordOperator(op0, n)
                s = random_stabilizer(n)
                @test apply!(copy(s),op)==apply!(copy(s),SingleQubitOperator(op))==apply!(copy(s),op_cc,[n])==apply!(copy(s),op_c)
                @test ==(apply!(copy(s),op,phases=false),apply!(copy(s),op_cc,[n],phases=false), phases=false)
                @test apply!(copy(s),op0)==apply!(copy(s),op0_c)
            end
            i = n÷2+1
            @test apply!(copy(s),sX(i)) == apply_single_x!(copy(s),i)
            @test apply!(copy(s),sY(i)) == apply_single_y!(copy(s),i)
            @test apply!(copy(s),sZ(i)) == apply_single_z!(copy(s),i)
            n==1 && continue
            s = random_stabilizer(n)
            i1,i2 = randperm(n)[1:2]
            @test apply!(copy(s),tCNOT,[i1,i2]) == apply!(copy(s),sCNOT(i1,i2))
            @test apply!(copy(s),tSWAP,[i1,i2]) == apply!(copy(s),sSWAP(i1,i2))
            @test apply!(copy(s),tCPHASE,[i1,i2]) == apply!(copy(s),sCPHASE(i1,i2))
        end
        @test_throws DimensionMismatch SingleQubitOperator(tCNOT,1)
        @test_throws DimensionMismatch CliffordOperator(sHadamard(5),2)
        @test_throws ArgumentError CliffordOperator(sHadamard(5),6,compact=true)
    end

    @testset "Convert between small ops" begin
        for op in subtypes(QuantumClifford.AbstractSingleQubitOperator)
            op == SingleQubitOperator && continue
            sop = op(1)
            sqsop = SingleQubitOperator(sop)
            cop = CliffordOperator(sop,1)
            csqop = CliffordOperator(sqsop,1)
            tcop = CliffordOperator(op)
            stcop = SingleQubitOperator(tcop)
            s = random_destabilizer(1)
            @test sop*s == sqsop*s == cop*s == csqop*s == tcop*s == stcop*s
        end
        for op in subtypes(QuantumClifford.AbstractTwoQubitOperator)
            sop = op(1,2)
            cop = CliffordOperator(sop,2)
            ccop = CliffordOperator(sop,2; compact=true)
            @test_throws DimensionMismatch CliffordOperator(op(1,4),3)
            @test_throws ArgumentError CliffordOperator(sop,3; compact=true)
            tcop = CliffordOperator(op)
            s = random_destabilizer(2)
            @test sop*s == cop*s == tcop*s == ccop*s
        end
        for op in subtypes(QuantumClifford.AbstractTwoQubitOperator)
            sop = op(1,10)
            op1 = CliffordOperator(sop,20)
            op2 = sSWAP(2,10)*(CliffordOperator(sop,2; compact=true)⊗tensor_pow(tId1, 18))*CliffordOperator(sSWAP(2,10),20)
            @test op1 == op2
        end
    end

    @testset "SingleQubitOperator inv methods" begin
        for gate_type in [sHadamard, sX, sY, sZ, sId1, sPhase, sInvPhase, sHadamardXY, sHadamardYZ, sV, sInvV, sSQRTY, sInvSQRTY]
            n = rand(1:10)
            @test CliffordOperator(inv(SingleQubitOperator(gate_type(n))), n) == inv(CliffordOperator(gate_type(n), n))
            @test CliffordOperator(inv(gate_type(n)), n) == inv(CliffordOperator(gate_type(n), n))
        end
        for i in 1:10
            random_op = random_clifford1(i)
            @test CliffordOperator(inv(random_op), i) == inv(CliffordOperator(random_op, i))
            @test CliffordOperator(inv(SingleQubitOperator(random_op)), i) == inv(CliffordOperator(random_op, i))
        end
    end
end
