# TODO split in separate files
@testitem "Noisy" begin
    using Random

    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    using QuantumClifford.Experimental.NoisyCircuits

    import AbstractAlgebra

    @testset "SparseGate" begin
        g = SparseGate(random_clifford(2), randperm(10)[1:2])
        gi = inv(g)
        c = random_stabilizer(10)
        @assert apply!(apply!(copy(c), g), gi) == c
    end

    @testset "Noisy Gates" begin
        g1 = SparseGate(tId1, [1])
        g2 = SparseGate(tCNOT, [2,3])
        g3 = sCNOT(4,5)
        g4 = sHadamard(6)
        n = UnbiasedUncorrelatedNoise(1)
        ng1 = NoisyGate(g1, n)
        ng2 = NoisyGate(g2, n)
        ng3 = NoisyGate(g3, n)
        ng4 = NoisyGate(g4, n)
        ng5 = NoiseOp(n,[7])
        state = ghz(7)
        res1, _ = mctrajectory!(copy(state), [ng1,ng2,ng3,ng4,ng5])
        res2, _ = mctrajectory!(copy(state), [g1,g2,g3,g4])
        @test res1 != res2 # has a very small chance of failing
        resp = petrajectories(copy(state), [ng1,ng2,ng3,ng4,ng5])
        @test all(values(resp).==0)
    end

    @testset "Monte Carlo Purification examples" begin
        g1 = SparseGate(tCNOT, [1,3])
        g2 = SparseGate(tCNOT, [2,4])
        m = BellMeasurement([sMX(3),sMX(4)])
        good_bell_state = S"XX
        ZZ"
        canonicalize_rref!(good_bell_state)
        v = VerifyOp(good_bell_state,[1,2])
        n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.03))
        init = Register(MixedDestabilizer(good_bell_state⊗good_bell_state))
        with_purification = mctrajectories(init, [n,g1,g2,m,v], trajectories=500)
        @test with_purification[failure_stat] > 5
        @test with_purification[false_success_stat] > 10
        @test with_purification[true_success_stat] > 420
        without_purification = mctrajectories(init, [n,v], trajectories=500)
        @test get(without_purification,failure_stat,0) == 0
        @test without_purification[false_success_stat] > 10
        @test without_purification[true_success_stat] > 450
        nonoise = mctrajectories(init, [g1,g2,m,v], trajectories=10)
        @test get(nonoise,failure_stat,0) == 0
        @test get(nonoise,false_success_stat,0) == 0
        @test nonoise[true_success_stat] == 10
    end

    @testset "Perturbative expansion Purification examples" begin
        @testset "Comparison to MC" begin
            compare(a,b, symbol) = abs(a[symbol]/500-b[symbol]) / (a[symbol]/500+b[symbol]+1e-5) < 0.3
            g1 = SparseGate(tCNOT, [1,3])
            g2 = SparseGate(tCNOT, [2,4])
            m = BellMeasurement([sMX(3),sMX(4)])
            good_bell_state = S"XX
            ZZ"
            canonicalize_rref!(good_bell_state)
            v = VerifyOp(good_bell_state,[1,2])
            n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.03))
            init = Register(MixedDestabilizer(good_bell_state⊗good_bell_state))
            mc = mctrajectories(init, [n,g1,g2,m,v], trajectories=500)
            pe = petrajectories(init, [n,g1,g2,m,v])
            @test compare(mc,pe,failure_stat)
            @test compare(mc,pe,false_success_stat)
            @test compare(mc,pe,true_success_stat)
            mc = mctrajectories(init, [n,v], trajectories=500)
            pe = petrajectories(init, [n,v])
            @test compare(mc,pe,failure_stat)
            @test compare(mc,pe,false_success_stat)
            @test compare(mc,pe,true_success_stat)
            mc = mctrajectories(init, [g1,g2,m,v], trajectories=500)
            pe = petrajectories(init, [g1,g2,m,v])
            @test compare(mc,pe,failure_stat)
            @test compare(mc,pe,false_success_stat)
            @test compare(mc,pe,true_success_stat)
        end

        @testset "Symbolic" begin
            R, (e,) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.RealField, ["e"])
            unity = R(1);
            good_bell_state = Register(MixedDestabilizer(S"XX ZZ"))
            initial_state = good_bell_state⊗good_bell_state
            g1 = SparseGate(tCNOT, [1,3]) # CNOT between qubit 1 and qubit 3 (both with Alice)
            g2 = SparseGate(tCNOT, [2,4]) # CNOT between qubit 2 and qubit 4 (both with Bob)
            m = BellMeasurement([sMX(3),sMX(4)]) # Bell measurement on qubit 3 and 4
            v = VerifyOp(good_bell_state,[1,2]) # Verify that qubit 1 and 2 indeed form a good Bell pair
            epsilon = e # The X or Y or Z error rate
            n = NoiseOpAll(UnbiasedUncorrelatedNoise(3epsilon))
            circuit = [n,g1,g2,m,v]
            pe_symbolic = petrajectories(initial_state, circuit, branch_weight=unity) # perturbative expansion
            @test pe_symbolic[false_success_stat] == -162.0*e^4 + 162.0*e^3 + -54.0*e^2 + 6.0*e
            @test pe_symbolic[failure_stat]   == -108.0*e^4 + 108.0*e^3 + -36.0*e^2 + 4.0*e
            @test pe_symbolic[true_success_stat]       == 27.0*e^4 + -54.0*e^3 + 36.0*e^2 + -10.0*e + 1.0
        end
    end

    @testset "Measurements" begin
        @testset "BellMeasurements" begin
            stateX = S"X"
            mX = BellMeasurement([sMX(1)])
            vX = VerifyOp(S"X", [1])
            determinate1 = mctrajectories(Register(MixedDestabilizer(stateX)), [mX,vX], trajectories=10)
            @test determinate1[failure_stat] == 0
            @test determinate1[false_success_stat] == 0
            @test determinate1[true_success_stat] == 10
            determinate1_pe = petrajectories(Register(MixedDestabilizer(copy(stateX))), [mX,vX])
            @test determinate1_pe[failure_stat] == 0
            @test determinate1_pe[false_success_stat] == 0
            @test determinate1_pe[true_success_stat] == 1
            stateZ = S"Z"
            random1 = mctrajectories(Register(MixedDestabilizer(stateZ)), [mX,vX], trajectories=500)
            @test random1[failure_stat] > 200
            @test random1[false_success_stat] == 0
            @test random1[true_success_stat] > 200
            random1_pe = petrajectories(Register(MixedDestabilizer(copy(stateZ))), [mX,vX])
            @test random1_pe[failure_stat] > 0.4
            @test random1_pe[false_success_stat] == 0
            @test random1_pe[true_success_stat] > 0.4
            bell_state = S" XX
            ZZ"
            m1 = BellMeasurement([sMX(1),sMX(2)])
            determinate2 = mctrajectories(Register(bell_state), [m1], trajectories=10)
            @test determinate2[failure_stat] == 0
            @test determinate2[false_success_stat] == 0
            @test determinate2[continue_stat] == 10
            determinate2_pe = petrajectories(Register(copy(bell_state)), [m1])
            @test determinate2_pe[failure_stat] == 0
            @test determinate2_pe[false_success_stat] == 0
            m2 = BellMeasurement([sMX(1),sMZ(2)])
            v = VerifyOp(bell_state, [1,2])
            random2 = mctrajectories(Register(bell_state), [m2,v], trajectories=500)
            @test random2[failure_stat]+random2[false_success_stat] == 500
            @test random2[true_success_stat] == 0
            random2_pe = petrajectories(Register(copy(bell_state)), [m2,v])
            @test random2_pe[failure_stat]+random2_pe[false_success_stat] == 1
            @test random2_pe[true_success_stat] == 0
        end

        @testset "PauliMeasurements" begin
            ghzState = S"XXX
            ZZI
            IZZ"
            m1 = PauliMeasurement(P"ZZI", 1)
            v = VerifyOp(ghzState, [1,2,3])
            register1 = Register(ghzState, zeros(Bool, 1))
            determinate1 = mctrajectories(register1, [m1,v], trajectories=10)
            @test determinate1[failure_stat] == 0
            @test determinate1[false_success_stat] == 0
            @test determinate1[true_success_stat] == 10
            determinate1_pe = petrajectories(register1, [m1,v])
            @test determinate1_pe[failure_stat] == 0
            @test determinate1_pe[false_success_stat] == 0
            @test determinate1_pe[true_success_stat] == 1
            m2 = PauliMeasurement(P"ZII", 1)
            register1 = Register(ghzState, zeros(Bool, 1))
            random1 = mctrajectories(register1, [m2,v], trajectories=50)
            @test random1[failure_stat] == 0
            @test random1[false_success_stat] == 50
            @test random1[true_success_stat] == 0
            random1_pe = petrajectories(register1, [m2,v])
            @test random1_pe[failure_stat] == 0
            @test random1_pe[false_success_stat] == 1
            @test random1_pe[true_success_stat] == 0
            m3 = PauliMeasurement(P"XII", 1)
            register1 = Register(ghzState, zeros(Bool, 1))
            random2 = mctrajectories(register1, [m3,v], trajectories=50)
            @test random2[failure_stat] == 0
            @test random2[false_success_stat] == 50
            @test random2[true_success_stat] == 0
            random2_pe = petrajectories(register1, [m3,v])
            @test random2_pe[failure_stat] == 0
            @test random2_pe[false_success_stat] == 1
            @test random2_pe[true_success_stat] == 0
        end

        @testset "Sparse Measurements" begin
            ghzState = S"XXX
            ZZI
            IZZ"
            v = VerifyOp(ghzState, [1,2,3])
            #= TODO reintroduce this type of SparseMeasurement (more than one qubit (so not sMZ), but not all qubits (so not DenseMeasurement))
            m1 = SparseMeasurement(P"ZZ", [1,2], 1)
            register1 = Register(ghzState, zeros(Bool, 1))
            determinate1 = mctrajectories(register1, [m1,v], trajectories=10)
            @test determinate1[failure_stat] == 0
            @test determinate1[false_success_stat] == 0
            @test determinate1[true_success_stat] == 10
            determinate1_pe = petrajectories(register1, [m1,v])
            @test determinate1_pe[failure_stat] == 0
            @test determinate1_pe[false_success_stat] == 0
            @test determinate1_pe[true_success_stat] == 1
            =#
            m2 = sMZ(1, 1)
            register1 = Register(ghzState, zeros(Bool, 1))
            random1 = mctrajectories(register1, [m2,v], trajectories=50)
            @test random1[failure_stat] == 0
            @test random1[false_success_stat] == 50
            @test random1[true_success_stat] == 0
            random1_pe = petrajectories(register1, [m2,v])
            @test random1_pe[failure_stat] == 0
            @test random1_pe[false_success_stat] == 1
            @test random1_pe[true_success_stat] == 0
            m3 = sMX(1, 1)
            register1 = Register(ghzState, zeros(Bool, 1))
            random2 = mctrajectories(register1, [m3,v], trajectories=50)
            @test random2[failure_stat] == 0
            @test random2[false_success_stat] == 50
            @test random2[true_success_stat] == 0
            random2_pe = petrajectories(register1, [m3,v])
            @test random2_pe[failure_stat] == 0
            @test random2_pe[false_success_stat] == 1
            @test random2_pe[true_success_stat] == 0
        end

        @testset "Conforming to the project! interface" begin
            state = Register(MixedDestabilizer(S"ZZ"), zeros(Bool, 1))
            meas = PauliMeasurement(P"ZI", 1)
            state, flag = applywstatus!(state, meas)
            @test rank(state.stab) == 2
            tab(state.stab).phases .= 0
            @test stabilizerview(state.stab) == S"ZZ
            ZI"
        end
    end

    @testset "Classical Bits" begin
        @testset "DecisionGate" begin
            X_error = CliffordOperator([P"X", P"-Z"])
            # testing single digit return value from decision function
            for s in [S"Z", S"-Z", S"X", S"-X", S"Y", S"-Y"]
                r = Register(s, [false])
                applywstatus!(r, PauliMeasurement(P"Z", 1))
                correctiveGate = SparseGate(X_error, [1])
                decisionFunction = syndrome -> syndrome[1] ? 1 : nothing
                applywstatus!(r, DecisionGate([correctiveGate], decisionFunction))
                @test stabilizerview(r) == S"Z"
            end

            # testing an array return from decision function
            expectedFinalState = S"ZI
            IZ"
            s = QuantumClifford.bell()
            r = Register(s, [false])
            applywstatus!(r, PauliMeasurement(P"ZI", 1))
            correctiveGates = [SparseGate(X_error, [1]), SparseGate(X_error, [2])]
            decisionFunction = syndrome -> syndrome[1] ? [1,2] : nothing
            applywstatus!(r, DecisionGate(correctiveGates, decisionFunction))
            canonicalize!(quantumstate(r))
            @test stabilizerview(r) == expectedFinalState

            s = QuantumClifford.bell((false, true)) # |01>+|10>
            r = Register(s, [false])
            applywstatus!(r, sMZ(1, 1))
            # we use the same corrective gates, with a different decision function
            decisionFunction = syndrome -> syndrome[1] ? [1] : 2 # both [1] and 1 should work
            applywstatus!(r, DecisionGate(correctiveGates, decisionFunction))
            canonicalize!(quantumstate(r))
            @test stabilizerview(r) == expectedFinalState
        end

        @testset "ConditionalGate" begin
            id_op = CliffordOperator([P"X", P"Z"])
            X_error = CliffordOperator([P"X", P"-Z"])

            for s in [S"Z", S"-Z", S"X", S"-X", S"Y", S"-Y"]
                r = Register(s, [false])
                applywstatus!(r, PauliMeasurement(P"Z", 1))
                correctiveGate = SparseGate(X_error, [1])
                identityGate = SparseGate(id_op, [1])
                applywstatus!(r, ConditionalGate(correctiveGate, identityGate, 1))
                @test stabilizerview(r) == S"Z"
            end

            expectedFinalState = S"ZI
            IZ"
            s = QuantumClifford.bell((false, true))
            r = Register(s, [false])
            applywstatus!(r, PauliMeasurement(P"ZI", 1))
            correctiveGate1 = SparseGate(X_error, [1])
            correctiveGate2 = SparseGate(X_error, [2])
            applywstatus!(r, ConditionalGate(correctiveGate1, correctiveGate2, 1))
            canonicalize!(quantumstate(r))
            @test stabilizerview(r) == expectedFinalState
        end
    end

    @testset "VerifyOp" begin
        using QuantumClifford.ECC: Steane7, parity_checks, naive_encoding_circuit
        @testset "Stabilizer passed as good_state is not a logical state" begin
            good_state = parity_checks(Steane7()) #passing in a code instead of a state within codespace
            verify = VerifyOp(good_state, 1:7)
            reg = Register(one(MixedDestabilizer,7),6) #dummy register to pass into applywstatus!
        
            
            @test_throws ArgumentError applywstatus!(reg, verify) #should throw an error since good_state isn't a logical state
            
    
        end
    
        @testset "Accepts pure good_state argument" begin
            good_state = S"ZZI
                            IZZ 
                            XXX" # passing in GHZ state as good_state argument since it exists within the repetition code's codespace
            verify = VerifyOp(good_state, 1:3)
            reg = Register(one(MixedDestabilizer,3),2) 
            encoding = [sHadamard(1)  ,  sCNOT(1,2)   ,sCNOT(1,3) ] # encoding the register's state into a logical state within the repetition code's codespace
            for i in encoding
                apply!(reg, i)
            end
    
            _, status = applywstatus!(reg, verify)
            @test status == true_success_stat
        end
    
      
    end
end
