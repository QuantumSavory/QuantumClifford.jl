using QuantumClifford
using Test

@testset "syndrome measurement of 3 qubit repetition code" begin
    circuit = [
        sX(1),
        sCNOT(1,2), sCNOT(1,3), # encode
        PauliError(2,0.75), # error
        sCNOT(1,4), sCNOT(2,4), sCNOT(2,5), sCNOT(3,5), sMZ(4,1), sMZ(5,2) # syndrome measurement
    ]
    frame = PauliFrame(100, 5, 2)
    frame = pftrajectories(frame, circuit)
    m = pfmeasurements(frame)
    # If the x component is set on the second qubit, then the fourth and fifth qubits should also have it set
    for (i, row) in enumerate(frame.frame)
        if row[2] == (true, false)
            @test row[4][1] && row[5][1] && m[i,1] && m[i,2]
        end
    end
end

@testset "GHZ correlations" begin
    ghz_circuit = [
        sHadamard(1), sCNOT(1,2), sCNOT(1,3), # prepare a GHZ state
        sMZ(1,1), sMZ(2,2), sMZ(3,3) # measure each qubit
    ]
    frame = PauliFrame(10^6, 3, 3)
    f = pftrajectories(frame, ghz_circuit)
    m = pfmeasurements(f)
    rowtotal_1s = sum(m, dims=2)[:,1]
    @test all(rowtotal_1s .% 3 .== 0)
    fractotal_1s = sum(rowtotal_1s)/3 / 10^6
    @test (fractotal_1s > 0.49) && (fractotal_1s < 0.51)
end

@testset "sMZ vs sMZR - mctrajectories vs pftrajectories" begin
    state = Register(one(MixedDestabilizer, 3), 6)
    frame = PauliFrame(100, 3, 6)

    ghz_circuit1 = [
        sHadamard(1), sCNOT(1,2), sCNOT(1,3), # prepare a GHZ state
        sMZ(1,1), sMZ(2,2), sMZ(3,3), # measure each qubit
        sMZ(1,4), sMZ(2,5), sMZ(3,6)  # measure each qubit again
    ]
    ms1 = stack([bitview(mctrajectory!(copy(state), ghz_circuit1)[1]) for i in 1:100], dims=1)
    mf1 = pfmeasurements(pftrajectories(copy(frame), ghz_circuit1))
    @test all(0.25 .< sum(ms1, dims=1)./100 .< 0.75)
    @test all(0.25 .< sum(mf1, dims=1)./100 .< 0.75)

    ghz_circuit2 = [
        sHadamard(1), sCNOT(1,2), sCNOT(1,3), # prepare a GHZ state
        sMRZ(1,1), sMZ(2,2), sMZ(3,3), # measure and reset each qubit
        sMZ(1,4), sMZ(2,5), sMZ(3,6)  # measure each qubit again
    ]
    ms2 = stack([bitview(mctrajectory!(copy(state), ghz_circuit2)[1]) for i in 1:100], dims=1)
    mf2 = pfmeasurements(pftrajectories(copy(frame), ghz_circuit2))
    @test all(0.25.*[1 1 1 0 1 1] .<= sum(ms2, dims=1)./100 .<= 0.75.*[1 1 1 0 1 1])
    @test all(0.25.*[1 1 1 0 1 1] .<= sum(mf2, dims=1)./100 .<= 0.75.*[1 1 1 0 1 1])
end
