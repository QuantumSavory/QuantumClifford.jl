# should merge into other test file

@testitem "Apply Inv" begin
    @testset "Apply Inv" begin
       
        @test apply_inv!(S"X", sHadamard(1)) == apply!(S"X", inv(CliffordOperator(sHadamard(1), 1)))
        @test apply_inv!(S"Y", sHadamard(1)) == apply!(S"Y", inv(CliffordOperator(sHadamard(1), 1)))
        @test apply_inv!(S"Z", sHadamard(1)) == apply!(S"Z", inv(CliffordOperator(sHadamard(1), 1)))

        @test apply_inv!(S"X", sHadamardXY(1)) == apply!(S"X", inv(CliffordOperator(sHadamardXY(1), 1)))
        @test apply_inv!(S"Y", sHadamardXY(1)) == apply!(S"Y", inv(CliffordOperator(sHadamardXY(1), 1)))
        @test apply_inv!(S"Z", sHadamardXY(1)) == apply!(S"Z", inv(CliffordOperator(sHadamardXY(1), 1)))

        @test apply_inv!(S"X", sHadamardYZ(1)) == apply!(S"X", inv(CliffordOperator(sHadamardYZ(1), 1)))
        @test apply_inv!(S"Y", sHadamardYZ(1)) == apply!(S"Y", inv(CliffordOperator(sHadamardYZ(1), 1)))
        @test apply_inv!(S"Z", sHadamardYZ(1)) == apply!(S"Z", inv(CliffordOperator(sHadamardYZ(1), 1)))

        @test apply_inv!(S"X", sPhase(1)) == apply!(S"X", inv(CliffordOperator(sPhase(1), 1)))
        @test apply_inv!(S"Y", sPhase(1)) == apply!(S"Y", inv(CliffordOperator(sPhase(1), 1)))
        @test apply_inv!(S"Z", sPhase(1)) == apply!(S"Z", inv(CliffordOperator(sPhase(1), 1)))

        @test apply_inv!(S"X", sInvPhase(1)) == apply!(S"X", inv(CliffordOperator(sInvPhase(1), 1)))
        @test apply_inv!(S"Y", sInvPhase(1)) == apply!(S"Y", inv(CliffordOperator(sInvPhase(1), 1)))
        @test apply_inv!(S"Z", sInvPhase(1)) == apply!(S"Z", inv(CliffordOperator(sInvPhase(1), 1)))

        @test apply_inv!(S"X", sX(1)) == apply!(S"X", inv(CliffordOperator(sX(1), 1)))
        @test apply_inv!(S"Y", sX(1)) == apply!(S"Y", inv(CliffordOperator(sX(1), 1)))
        @test apply_inv!(S"Z", sX(1)) == apply!(S"Z", inv(CliffordOperator(sX(1), 1)))

        @test apply_inv!(S"X", sY(1)) == apply!(S"X", inv(CliffordOperator(sY(1), 1)))
        @test apply_inv!(S"Y", sY(1)) == apply!(S"Y", inv(CliffordOperator(sY(1), 1)))
        @test apply_inv!(S"Z", sY(1)) == apply!(S"Z", inv(CliffordOperator(sY(1), 1)))

        @test apply_inv!(S"X", sZ(1)) == apply!(S"X", inv(CliffordOperator(sZ(1), 1)))
        @test apply_inv!(S"Y", sZ(1)) == apply!(S"Y", inv(CliffordOperator(sZ(1), 1)))
        @test apply_inv!(S"Z", sZ(1)) == apply!(S"Z", inv(CliffordOperator(sZ(1), 1)))

        @test apply_inv!(S"X", sSQRTX(1)) == apply!(S"X", inv(CliffordOperator(sSQRTX(1), 1)))
        @test apply_inv!(S"Y", sSQRTX(1)) == apply!(S"Y", inv(CliffordOperator(sSQRTX(1), 1)))
        @test apply_inv!(S"Z", sSQRTX(1)) == apply!(S"Z", inv(CliffordOperator(sSQRTX(1), 1)))

        @test apply_inv!(S"X", sInvSQRTX(1)) == apply!(S"X", inv(CliffordOperator(sInvSQRTX(1), 1)))
        @test apply_inv!(S"Y", sInvSQRTX(1)) == apply!(S"Y", inv(CliffordOperator(sInvSQRTX(1), 1)))
        @test apply_inv!(S"Z", sInvSQRTX(1)) == apply!(S"Z", inv(CliffordOperator(sInvSQRTX(1), 1)))

        @test apply_inv!(S"X", sSQRTY(1)) == apply!(S"X", inv(CliffordOperator(sSQRTY(1), 1)))
        @test apply_inv!(S"Y", sSQRTY(1)) == apply!(S"Y", inv(CliffordOperator(sSQRTY(1), 1)))
        @test apply_inv!(S"Z", sSQRTY(1)) == apply!(S"Z", inv(CliffordOperator(sSQRTY(1), 1)))

        @test apply_inv!(S"X", sInvSQRTY(1)) == apply!(S"X", inv(CliffordOperator(sInvSQRTY(1), 1)))
        @test apply_inv!(S"Y", sInvSQRTY(1)) == apply!(S"Y", inv(CliffordOperator(sInvSQRTY(1), 1)))
        @test apply_inv!(S"Z", sInvSQRTY(1)) == apply!(S"Z", inv(CliffordOperator(sInvSQRTY(1), 1)))

        @test apply_inv!(S"X", sCXYZ(1)) == apply!(S"X", inv(CliffordOperator(sCXYZ(1), 1)))
        @test apply_inv!(S"Y", sCXYZ(1)) == apply!(S"Y", inv(CliffordOperator(sCXYZ(1), 1)))
        @test apply_inv!(S"Z", sCXYZ(1)) == apply!(S"Z", inv(CliffordOperator(sCXYZ(1), 1)))

        @test apply_inv!(S"X", sCZYX(1)) == apply!(S"X", inv(CliffordOperator(sCZYX(1), 1)))
        @test apply_inv!(S"Y", sCZYX(1)) == apply!(S"Y", inv(CliffordOperator(sCZYX(1), 1)))
        @test apply_inv!(S"Z", sCZYX(1)) == apply!(S"Z", inv(CliffordOperator(sCZYX(1), 1)))

    end
end