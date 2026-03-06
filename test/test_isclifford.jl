@testitem "isclifford" begin
    using QuantumClifford

    @testset "Single-qubit symbolic Cliffords" begin
        @test isclifford(sHadamard(1)) == true
        @test isclifford(sPhase(1)) == true
        @test isclifford(sInvPhase(1)) == true
        @test isclifford(sX(1)) == true
        @test isclifford(sY(1)) == true
        @test isclifford(sZ(1)) == true
        @test isclifford(sId1(1)) == true
        @test isclifford(sSQRTX(1)) == true
        @test isclifford(sInvSQRTX(1)) == true
        @test isclifford(sSQRTY(1)) == true
        @test isclifford(sInvSQRTY(1)) == true
        @test isclifford(sHadamardXY(1)) == true
        @test isclifford(sHadamardYZ(1)) == true
        @test isclifford(sCXYZ(1)) == true
        @test isclifford(sCZYX(1)) == true
    end

    @testset "Two-qubit symbolic Cliffords" begin
        @test isclifford(sCNOT(1, 2)) == true
        @test isclifford(sCPHASE(1, 2)) == true
        @test isclifford(sSWAP(1, 2)) == true
        @test isclifford(sXCX(1, 2)) == true
        @test isclifford(sXCY(1, 2)) == true
        @test isclifford(sXCZ(1, 2)) == true
        @test isclifford(sYCX(1, 2)) == true
        @test isclifford(sYCY(1, 2)) == true
        @test isclifford(sYCZ(1, 2)) == true
        @test isclifford(sZCX(1, 2)) == true
        @test isclifford(sZCY(1, 2)) == true
        @test isclifford(sZCZ(1, 2)) == true
        @test isclifford(sSWAPCX(1, 2)) == true
        @test isclifford(sInvSWAPCX(1, 2)) == true
        @test isclifford(sCZSWAP(1, 2)) == true
        @test isclifford(sCXSWAP(1, 2)) == true
        @test isclifford(sISWAP(1, 2)) == true
        @test isclifford(sInvISWAP(1, 2)) == true
        @test isclifford(sSQRTZZ(1, 2)) == true
        @test isclifford(sInvSQRTZZ(1, 2)) == true
    end

    @testset "Measurements" begin
        @test isclifford(sMX(1, 1)) == true
        @test isclifford(sMY(1, 1)) == true
        @test isclifford(sMZ(1, 1)) == true
        @test isclifford(sMRX(1, 1)) == true
        @test isclifford(sMRY(1, 1)) == true
        @test isclifford(sMRZ(1, 1)) == true
        @test isclifford(PauliMeasurement(P"X", 1)) == true
    end

    @testset "SparseGate" begin
        sparse_clifford = SparseGate(CliffordOperator(tHadamard), [1])
        @test isclifford(sparse_clifford) == true
    end

    @testset "Dense CliffordOperator" begin
        @test isclifford(CliffordOperator(tHadamard)) == true
    end
end
