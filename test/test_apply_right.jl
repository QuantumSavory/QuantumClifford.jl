@testitem "Apply Right" begin
    using InteractiveUtils
    using QuantumClifford: AbstractCliffordOperator
    
    # SLOW version of apply_right! for testing
    function apply_right_slow!(l::CliffordOperator, r::AbstractCliffordOperator; phases=true)
        apply!(CliffordOperator(r, nqubits(l)), l; phases=phases)
    end
    
    q = 64
    shots = 16

    # @testset "Apply Right single-qubit" begin
    #     for gate in subtypes(AbstractSingleQubitOperator)
    #         if gate in [sPhase, sInvPhase, SingleQubitOperator, sCXYZ, sCZYX]
    #             continue
    #         end
    #         r = gate(rand(1:q))
            
    #         for _ in 1:shots
    #             l = random_clifford(q)
    #             @test isequal(apply_right!(copy(l), r), apply_right_slow!(l, r))
    #         end
    #     end
    # end

    # @testset "Apply Right two-qubit" begin
    #     for gate in subtypes(AbstractTwoQubitOperator)
    #         if gate in [sSWAPCX, sInvSWAPCX, sCZSWAP, sCXSWAP, sCPHASE, sZCrY, sInvZCrY]
    #             continue
    #         end
    #         q1 = rand(1:q); q2 = rand(setdiff(1:q, [q1]))
    #         r = gate(q1, q2)

    #         for _ in 1:shots
    #             l = random_clifford(q)
    #             @test isequal(apply_right!(copy(l), r), apply_right_slow!(l, r))
    #         end
    #     end
    # end

    @testset "Apply Right single-qubit" begin
        for _ in 1:shots
            l = random_clifford(q)
            q1 = rand(1:q)

            @test isequal(apply_right!(copy(l), sHadamard(q1)), apply_right_slow!(l, sHadamard(q1)))
            @test isequal(apply_right!(copy(l), sHadamardXY(q1)), apply_right_slow!(l, sHadamardXY(q1)))
            @test isequal(apply_right!(copy(l), sHadamardYZ(q1)), apply_right_slow!(l, sHadamardYZ(q1)))
            # @test isequal(apply_right!(copy(l), sPhase(q1)), apply_right_slow!(l, sPhase(q1)))
            # @test isequal(apply_right!(copy(l), sInvPhase(q1)), apply_right_slow!(l, sInvPhase(q1)))
            @test isequal(apply_right!(copy(l), sX(q1)), apply_right_slow!(l, sX(q1)))
            @test isequal(apply_right!(copy(l), sY(q1)), apply_right_slow!(l, sY(q1)))
            @test isequal(apply_right!(copy(l), sZ(q1)), apply_right_slow!(l, sZ(q1)))
            @test isequal(apply_right!(copy(l), sSQRTX(q1)), apply_right_slow!(l, sSQRTX(q1)))
            @test isequal(apply_right!(copy(l), sInvSQRTX(q1)), apply_right_slow!(l, sInvSQRTX(q1)))
            @test isequal(apply_right!(copy(l), sSQRTY(q1)), apply_right_slow!(l, sSQRTY(q1)))
            @test isequal(apply_right!(copy(l), sInvSQRTY(q1)), apply_right_slow!(l, sInvSQRTY(q1)))
            # @test isequal(apply_right!(copy(l), sSQRTZ(q1)), apply_right_slow!(l, sSQRTZ(q1)))
            # @test isequal(apply_right!(copy(l), sInvSQRTZ(q1)), apply_right_slow!(l, sInvSQRTZ(q1)))
            # @test isequal(apply_right!(copy(l), sCXYZ(q1)), apply_right_slow!(l, sCXYZ(q1)))
            # @test isequal(apply_right!(copy(l), sCZYX(q1)), apply_right_slow!(l, sCZYX(q1)))
            @test isequal(apply_right!(copy(l), sId1(q1)), apply_right_slow!(l, sId1(q1)))
        end
    end

    @testset "Apply Right two-qubit" begin
        for _ in 1:shots
            l = random_clifford(q)
            q1 = rand(1:q); q2 = rand(setdiff(1:q, [q1]))

            @test isequal(apply_right!(copy(l), sSWAP(q1, q2)), apply_right_slow!(l, sSWAP(q1, q2)))
            # @test isequal(apply_right!(copy(l), sSWAPCX(q1, q2)), apply_right_slow!(l, sSWAPCX(q1, q2)))
            # @test isequal(apply_right!(copy(l), sInvSWAPCX(q1, q2)), apply_right_slow!(l, sInvSWAPCX(q1, q2)))
            @test isequal(apply_right!(copy(l), sISWAP(q1, q2)), apply_right_slow!(l, sISWAP(q1, q2)))
            @test isequal(apply_right!(copy(l), sInvISWAP(q1, q2)), apply_right_slow!(l, sInvISWAP(q1, q2)))
            # @test isequal(apply_right!(copy(l), sCZSWAP(q1, q2)), apply_right_slow!(l, sCZSWAP(q1, q2)))
            # @test isequal(apply_right!(copy(l), sCXSWAP(q1, q2)), apply_right_slow!(l, sCXSWAP(q1, q2)))
            @test isequal(apply_right!(copy(l), sCNOT(q1, q2)), apply_right_slow!(l, sCNOT(q1, q2)))
            # @test isequal(apply_right!(copy(l), sCPHASE(q1, q2)), apply_right_slow!(l, sCPHASE(q1, q2)))
            @test isequal(apply_right!(copy(l), sZCX(q1, q2)), apply_right_slow!(l, sZCX(q1, q2)))
            @test isequal(apply_right!(copy(l), sZCY(q1, q2)), apply_right_slow!(l, sZCY(q1, q2)))
            @test isequal(apply_right!(copy(l), sZCZ(q1, q2)), apply_right_slow!(l, sZCZ(q1, q2)))
            @test isequal(apply_right!(copy(l), sXCX(q1, q2)), apply_right_slow!(l, sXCX(q1, q2)))
            @test isequal(apply_right!(copy(l), sXCY(q1, q2)), apply_right_slow!(l, sXCY(q1, q2)))
            @test isequal(apply_right!(copy(l), sXCZ(q1, q2)), apply_right_slow!(l, sXCZ(q1, q2)))
            @test isequal(apply_right!(copy(l), sYCX(q1, q2)), apply_right_slow!(l, sYCX(q1, q2)))
            @test isequal(apply_right!(copy(l), sYCY(q1, q2)), apply_right_slow!(l, sYCY(q1, q2)))
            @test isequal(apply_right!(copy(l), sYCZ(q1, q2)), apply_right_slow!(l, sYCZ(q1, q2)))
            # @test isequal(apply_right!(copy(l), sZCrY(q1, q2)), apply_right_slow!(l, sZCrY(q1, q2)))
            # @test isequal(apply_right!(copy(l), sInvZCrY(q1, q2)), apply_right_slow!(l, sInvZCrY(q1, q2)))
            @test isequal(apply_right!(copy(l), sSQRTZZ(q1, q2)), apply_right_slow!(l, sSQRTZZ(q1, q2)))
            @test isequal(apply_right!(copy(l), sInvSQRTZZ(q1, q2)), apply_right_slow!(l, sInvSQRTZZ(q1, q2)))
            @test isequal(apply_right!(copy(l), sSQRTXX(q1, q2)), apply_right_slow!(l, sSQRTXX(q1, q2)))
            @test isequal(apply_right!(copy(l), sInvSQRTXX(q1, q2)), apply_right_slow!(l, sInvSQRTXX(q1, q2)))
            @test isequal(apply_right!(copy(l), sSQRTYY(q1, q2)), apply_right_slow!(l, sSQRTYY(q1, q2)))
            @test isequal(apply_right!(copy(l), sInvSQRTYY(q1, q2)), apply_right_slow!(l, sInvSQRTYY(q1, q2)))
        end
    end
end