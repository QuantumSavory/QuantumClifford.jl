@testitem "Pauli frame" begin
    @testset "syndrome measurement of 3 qubit repetition code" begin
        circuit = [
                   sX(1),
                   sCNOT(1,2), sCNOT(1,3), # encode
                   PauliError(2,0.75), # error
                   sCNOT(1,4), sCNOT(2,4), sCNOT(2,5), sCNOT(3,5), sMZ(4,1), sMZ(5,2) # syndrome measurement
                  ]
        frame = PauliFrame(100, 5, 2)
        frame = pftrajectories(frame, circuit)
        frame1 = pftrajectories(circuit; trajectories=100, threads=false)
        frame2 = pftrajectories(circuit; trajectories=100, threads=true)
        # If the x component is set on the second qubit, then the fourth and fifth qubits should also have it set
        for f in [frame, frame1, frame2]
            m = pfmeasurements(f)
            for (i, row) in enumerate(f.frame)
                if row[2] == (true, false)
                    @test row[4][1] && row[5][1] && m[i,1] && m[i,2]
                end
            end
        end
    end

    @testset "GHZ correlations" begin
        ghz_circuit = [
                       sHadamard(1), sCNOT(1,2), sCNOT(1,3), # prepare a GHZ state
                       sMZ(1,1), sMZ(2,2), sMZ(3,3) # measure each qubit
                      ]
        n = 10^6
        frame = PauliFrame(n, 3, 3)
        f = pftrajectories(frame, ghz_circuit)
        m = pfmeasurements(f)
        frame1 = pftrajectories(ghz_circuit; trajectories=n, threads=false)
        m1 = pfmeasurements(frame1)
        frame2 = pftrajectories(ghz_circuit; trajectories=n, threads=true)
        m2 = pfmeasurements(frame2)
        for _m in [m, m1, m2]
            rowtotal_1s = sum(m, dims=2)[:,1]
            @test all(rowtotal_1s .% 3 .== 0)
            fractotal_1s = sum(rowtotal_1s)/3 / n
            @test (fractotal_1s > 0.49) && (fractotal_1s < 0.51)
        end
    end

    @testset "sMZ vs sMRZ - mctrajectories vs pftrajectories" begin
        n = 2000
        state = Register(one(MixedDestabilizer, 3), 6)
        frame = PauliFrame(n, 3, 6)

        ghz_circuit1 = [
                        sHadamard(1), sCNOT(1,2), sCNOT(1,3), # prepare a GHZ state
                        sMZ(1,1), sMZ(2,2), sMZ(3,3), # measure each qubit
                        sMZ(1,4), sMZ(2,5), sMZ(3,6)  # measure each qubit again
                       ]
        for m in [
                  stack([bitview(mctrajectory!(copy(state), ghz_circuit1)[1]) for i in 1:n], dims=1),
                  pfmeasurements(pftrajectories(copy(frame), ghz_circuit1)),
                  pfmeasurements(pftrajectories(ghz_circuit1;trajectories=n,threads=false)),
                  pfmeasurements(pftrajectories(ghz_circuit1;trajectories=n,threads=true)),
                 ]
            @test all(0.25 .< sum(m, dims=1)./n .< 0.75)
        end

        ghz_circuit2 = [
                        sHadamard(1), sCNOT(1,2), sCNOT(1,3), # prepare a GHZ state
                        sMRZ(1,1), sMZ(2,2), sMZ(3,3), # measure and reset each qubit
                        sMZ(1,4), sMZ(2,5), sMZ(3,6)  # measure each qubit again
                       ]
        for m in [
                  stack([bitview(mctrajectory!(copy(state), ghz_circuit2)[1]) for i in 1:n], dims=1),
                  pfmeasurements(pftrajectories(copy(frame), ghz_circuit2)),
                  pfmeasurements(pftrajectories(ghz_circuit2;trajectories=n,threads=false)),
                  pfmeasurements(pftrajectories(ghz_circuit2;trajectories=n,threads=true)),
                 ]
            @test all(0.25.*[1 1 1 0 1 1] .<= sum(m, dims=1)./n .<= 0.75.*[1 1 1 0 1 1])
        end

        noncom_circuit = [
                          sHadamard(1), sMRZ(1,1), sX(1), sMZ(1,2), sMRZ(1,3), sMRZ(1,4), sHadamard(1), sMZ(1,5)
                         ]
        ms3 = stack([bitview(mctrajectory!(copy(state), noncom_circuit)[1]) for i in 1:n], dims=1)
        @test all(0.25.*[1 4 4 0 1 0] .<= sum(ms3, dims=1)./n .<= 0.75.*[1 2 2 0 1 0])
        for m in [
                  pfmeasurements(pftrajectories(copy(frame), noncom_circuit)),
                  pfmeasurements(pftrajectories(noncom_circuit;trajectories=n,threads=false)),
                  pfmeasurements(pftrajectories(noncom_circuit;trajectories=n,threads=true)),
                 ]
            @test all(0.25.*[1 0 0 0 1] .<= (sum(m, dims=1)[:,1:5])./n .<= 0.75.*[1 0 0 0 1])
        end
    end

    @testset "PauliMeasurements" begin
        n = 2000
        state = Register(one(MixedDestabilizer, 3), 5)
        frame = PauliFrame(n, 3, 5)

        glassy_ghz_circuit = [
            sHadamard(1), sHadamard(2), sHadamard(3),
            PauliMeasurement(P"ZZ_", 1), PauliMeasurement(P"_ZZ", 2),
            sMZ(1, 3), sMZ(2, 4), sMZ(3, 5)
        ]
        for m in [pfmeasurements(pftrajectories(copy(frame), glassy_ghz_circuit)),
            pfmeasurements(pftrajectories(glassy_ghz_circuit; trajectories=n, threads=false)),
            pfmeasurements(pftrajectories(glassy_ghz_circuit; trajectories=n, threads=true))]

            # decode based on measurement outcomes
            for r in eachrow(m)
                r[4] ⊻= r[1]
                r[5] ⊻= r[1] ⊻ r[2]
            end

            # check that the correct correlations are present
            @test all(m[:, 3] .== m[:, 4] .== m[:, 5])
        end
    end
end
