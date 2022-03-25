using QuantumClifford.Experimental.NoisyCircuits
import AbstractAlgebra

function test_noisycircuits()
    @testset "Noisy Circuits" begin
        @testset "Monte Carlo sims" begin
            @testset "Purification examples" begin
                g1 = SparseGate(CNOT, [1,3])
                g2 = SparseGate(CNOT, [2,4])
                m = BellMeasurement([X,X],[3,4])
                good_bell_state = S"XX
                                    ZZ"
                canonicalize_rref!(good_bell_state)
                v = VerifyOp(good_bell_state,[1,2])
                n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.01))
                with_purification = mctrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v], trajectories=500)
                @test with_purification[:detected_failure] > 5
                @test with_purification[:undetected_failure] > 10
                @test with_purification[:true_success] > 430
                without_purification = mctrajectories(good_bell_state⊗good_bell_state, [n,v], trajectories=500)
                @test without_purification[:detected_failure] == 0
                @test without_purification[:undetected_failure] > 10
                @test without_purification[:true_success] > 450
                nonoise = mctrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v], trajectories=10)
                @test nonoise[:detected_failure] == 0
                @test nonoise[:undetected_failure] == 0
                @test nonoise[:true_success] == 10
            end
        end
        @testset "Perturbative expansion sims" begin
            @testset "Purification examples comparison to MC" begin
                compare(a,b, symbol) = abs(a[symbol]/500-b[symbol]) / (a[symbol]/500+b[symbol]+1e-5) < 0.3
                g1 = SparseGate(CNOT, [1,3])
                g2 = SparseGate(CNOT, [2,4])
                m = BellMeasurement([X,X],[3,4])
                good_bell_state = S"XX
                                    ZZ"
                canonicalize_rref!(good_bell_state)
                v = VerifyOp(good_bell_state,[1,2])
                n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.01))
                mc = mctrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v], trajectories=500)
                pe = petrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v])
                @test compare(mc,pe,:detected_failure)
                @test compare(mc,pe,:undetected_failure)
                @test compare(mc,pe,:true_success)
                mc = mctrajectories(good_bell_state⊗good_bell_state, [n,v], trajectories=500)
                pe = petrajectories(good_bell_state⊗good_bell_state, [n,v])
                @test compare(mc,pe,:detected_failure)
                @test compare(mc,pe,:undetected_failure)
                @test compare(mc,pe,:true_success)
                mc = mctrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v], trajectories=500)
                pe = petrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v])
                @test compare(mc,pe,:detected_failure)
                @test compare(mc,pe,:undetected_failure)
                @test compare(mc,pe,:true_success)
            end
            
            @testset "Symbolic" begin
                for statetype in [Stabilizer, MixedDestabilizer]
                    R, (e,) = AbstractAlgebra.PolynomialRing(AbstractAlgebra.RealField, ["e"])
                    unity = R(1);
                    
                    good_bell_state = statetype(S"XX
                                                ZZ")
                    initial_state = good_bell_state⊗good_bell_state

                    g1 = SparseGate(CNOT, [1,3]) # CNOT between qubit 1 and qubit 3 (both with Alice)
                    g2 = SparseGate(CNOT, [2,4]) # CNOT between qubit 2 and qubit 4 (both with Bob)
                    m = BellMeasurement([X,X],[3,4]) # Bell measurement on qubit 3 and 4
                    v = VerifyOp(good_bell_state,[1,2]) # Verify that qubit 1 and 2 indeed form a good Bell pair
                    epsilon = e # The error rate
                    n = NoiseOpAll(UnbiasedUncorrelatedNoise(epsilon))

                    # This circuit performs a depolarization at rate `epsilon` to all qubits,
                    # then bilater CNOT operations
                    # then a Bell measurement
                    # followed by checking whether the final result indeed corresponds to the correct Bell pair.
                    circuit = [n,g1,g2,m,v]

                    pe_symbolic = petrajectories(initial_state, circuit, branch_weight=unity) # perturbative expansion
                    @test pe_symbolic[:undetected_failure] == -162.0*e^4 + 162.0*e^3 + -54.0*e^2 + 6.0*e
                    @test pe_symbolic[:detected_failure]   == -108.0*e^4 + 108.0*e^3 + -36.0*e^2 + 4.0*e
                    @test pe_symbolic[:true_success]       == 27.0*e^4 + -54.0*e^3 + 36.0*e^2 + -10.0*e + 1.0
                end
            end
        end
        @testset "Measurements" begin
            # compare(a, b, symbols, c) = [abs(a[symbol]/c-b[symbol]) / (a[symbol]/c+b[symbol]+1e-5) < 0.1 for symbol in symbols]
            @testset "BellMeasurements" begin
                stateX = S"X"
                mX = BellMeasurement([X], [1])
                vX = VerifyOp(S"X", [1])
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    determinate1 = mctrajectories(stabType(stateX), [mX,vX], trajectories=10)
                    @test determinate1[:detected_failure] == 0
                    @test determinate1[:undetected_failure] == 0
                    @test determinate1[:true_success] == 10
                    determinate1_pe = petrajectories(stabType(copy(stateX)), [mX,vX])
                    @test determinate1_pe[:detected_failure] == 0
                    @test determinate1_pe[:undetected_failure] == 0
                    @test determinate1_pe[:true_success] == 1
                end
                stateZ = S"Z"
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    random1 = mctrajectories(stabType(stateZ), [mX,vX], trajectories=500)
                    @test random1[:detected_failure] > 200
                    @test random1[:undetected_failure] == 0
                    @test random1[:true_success] > 200
                    random1_pe = petrajectories(stabType(copy(stateZ)), [mX,vX])
                    @test random1_pe[:detected_failure] > 0.4
                    @test random1_pe[:undetected_failure] == 0
                    @test random1_pe[:true_success] > 0.4
                end
                bell_state = S" XX
                                ZZ"
                m1 = BellMeasurement([X,X], [1,2])
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    determinate2 = mctrajectories(stabType(bell_state), [m1], trajectories=10)
                    @test determinate2[:detected_failure] == 0
                    @test determinate2[:undetected_failure] == 0
                    @test determinate2[:continue] == 10
                    determinate2_pe = petrajectories(stabType(copy(bell_state)), [m1])
                    @test determinate2_pe[:detected_failure] == 0
                    @test determinate2_pe[:undetected_failure] == 0
                end
                m2 = BellMeasurement([X,Z], [1,2])
                v = VerifyOp(bell_state, [1,2])
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    random2 = mctrajectories(stabType(bell_state), [m2,v], trajectories=500)
                    @test random2[:detected_failure]+random2[:undetected_failure] == 500
                    @test random2[:true_success] == 0
                    random2_pe = petrajectories(stabType(copy(bell_state)), [m2,v])
                    @test random2_pe[:detected_failure]+random2_pe[:undetected_failure] == 1
                    @test random2_pe[:true_success] == 0
                end
            end
            @testset "DenseMeasurements" begin
                ghzState = S"XXX
                            ZZI
                            IZZ"
                m1 = DenseMeasurement(P"ZZI", 1)
                v = VerifyOp(ghzState, [1,2,3])
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    register1 = Register(stabType(ghzState), zeros(Bool, 1))
                    determinate1 = mctrajectories(register1, [m1,v], trajectories=10)
                    @test determinate1[:detected_failure] == 0
                    @test determinate1[:undetected_failure] == 0
                    @test determinate1[:true_success] == 10
                    determinate1_pe = petrajectories(register1, [m1,v])
                    @test determinate1_pe[:detected_failure] == 0
                    @test determinate1_pe[:undetected_failure] == 0
                    @test determinate1_pe[:true_success] == 1
                end
                m2 = DenseMeasurement(P"ZII", 1)
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    register1 = Register(stabType(ghzState), zeros(Bool, 1))
                    random1 = mctrajectories(register1, [m2,v], trajectories=50)
                    @test random1[:detected_failure] == 0
                    @test random1[:undetected_failure] == 50
                    @test random1[:true_success] == 0
                    random1_pe = petrajectories(register1, [m2,v])
                    @test random1_pe[:detected_failure] == 0
                    @test random1_pe[:undetected_failure] == 1
                    @test random1_pe[:true_success] == 0
                end
                m3 = DenseMeasurement(P"XII", 1)
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    register1 = Register(stabType(ghzState), zeros(Bool, 1))
                    random2 = mctrajectories(register1, [m3,v], trajectories=50)
                    @test random2[:detected_failure] == 0
                    @test random2[:undetected_failure] == 50
                    @test random2[:true_success] == 0
                    random2_pe = petrajectories(register1, [m3,v])
                    @test random2_pe[:detected_failure] == 0
                    @test random2_pe[:undetected_failure] == 1
                    @test random2_pe[:true_success] == 0
                end
            end
            @testset "SparseMeasurements" begin
                ghzState = S"XXX
                            ZZI
                            IZZ"
                m1 = SparseMeasurement(P"ZZ", [1,2], 1)
                v = VerifyOp(ghzState, [1,2,3])
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    register1 = Register(stabType(ghzState), zeros(Bool, 1))
                    determinate1 = mctrajectories(register1, [m1,v], trajectories=10)
                    @test determinate1[:detected_failure] == 0
                    @test determinate1[:undetected_failure] == 0
                    @test determinate1[:true_success] == 10
                    determinate1_pe = petrajectories(register1, [m1,v])
                    @test determinate1_pe[:detected_failure] == 0
                    @test determinate1_pe[:undetected_failure] == 0
                    @test determinate1_pe[:true_success] == 1
                end
                m2 = SparseMeasurement(P"Z", [1], 1)
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    register1 = Register(stabType(ghzState), zeros(Bool, 1))
                    random1 = mctrajectories(register1, [m2,v], trajectories=50)
                    @test random1[:detected_failure] == 0
                    @test random1[:undetected_failure] == 50
                    @test random1[:true_success] == 0
                    random1_pe = petrajectories(register1, [m2,v])
                    @test random1_pe[:detected_failure] == 0
                    @test random1_pe[:undetected_failure] == 1
                    @test random1_pe[:true_success] == 0
                end
                m3 = SparseMeasurement(P"X", [1], 1)
                for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
                    register1 = Register(stabType(ghzState), zeros(Bool, 1))
                    random2 = mctrajectories(register1, [m3,v], trajectories=50)
                    @test random2[:detected_failure] == 0
                    @test random2[:undetected_failure] == 50
                    @test random2[:true_success] == 0
                    random2_pe = petrajectories(register1, [m3,v])
                    @test random2_pe[:detected_failure] == 0
                    @test random2_pe[:undetected_failure] == 1
                    @test random2_pe[:true_success] == 0
                end
            end
            @testset "Conforming to the project! interface" begin
                state = Register(MixedDestabilizer(S"ZZ"), zeros(Bool, 1))
                meas = DenseMeasurement(P"ZI", 1)
                state, flag = applyop!(state, meas)
                @test state.stab.rank == 2
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
                    applyop!(r, DenseMeasurement(P"Z", 1))
                    correctiveGate = SparseGate(X_error, [1])
                    decisionFunction = syndrome -> syndrome[1] ? 1 : nothing
                    applyop!(r, DecisionGate([correctiveGate], decisionFunction))
                    @test stabilizerview(r) == S"Z"
                end

                # testing an array return from decision function
                expectedFinalState = S"ZI
                                    IZ"
                s = QuantumClifford.bell()
                r = Register(s, [false])
                applyop!(r, DenseMeasurement(P"ZI", 1))
                # applyop!(r, SparseMeasurement(P"Z", [1], 1))
                correctiveGates = [SparseGate(X_error, [1]), SparseGate(X_error, [2])]
                decisionFunction = syndrome -> syndrome[1] ? [1,2] : nothing
                applyop!(r, DecisionGate(correctiveGates, decisionFunction))
                canonicalize!(quantumstate(r))
                @test stabilizerview(r) == expectedFinalState

                s = QuantumClifford.bell((false, true)) # |01>+|10>
                r = Register(s, [false])
                # applyop!(r, DenseMeasurement(P"ZI", 1))
                applyop!(r, SparseMeasurement(P"Z", [1], 1))
                # we use the same corrective gates, with a different decision function
                decisionFunction = syndrome -> syndrome[1] ? [1] : 2 # both [1] and 1 should work
                applyop!(r, DecisionGate(correctiveGates, decisionFunction))
                canonicalize!(quantumstate(r))
                @test stabilizerview(r) == expectedFinalState
            end
            @testset "ConditionalGate" begin
                id_op = CliffordOperator([P"X", P"Z"])
                X_error = CliffordOperator([P"X", P"-Z"])

                for s in [S"Z", S"-Z", S"X", S"-X", S"Y", S"-Y"]
                    r = Register(s, [false])
                    applyop!(r, DenseMeasurement(P"Z", 1))
                    correctiveGate = SparseGate(X_error, [1])
                    identityGate = SparseGate(id_op, [1])
                    applyop!(r, ConditionalGate(correctiveGate, identityGate, 1))
                    @test stabilizerview(r) == S"Z"
                end

                expectedFinalState = S"ZI
                                       IZ"
                s = QuantumClifford.bell((false, true))
                r = Register(s, [false])
                applyop!(r, DenseMeasurement(P"ZI", 1))
                correctiveGate1 = SparseGate(X_error, [1])
                correctiveGate2 = SparseGate(X_error, [2])
                applyop!(r, ConditionalGate(correctiveGate1, correctiveGate2, 1))
                canonicalize!(quantumstate(r))
                @test stabilizerview(r) == expectedFinalState
            end
        end
    end
end

test_noisycircuits()