@testitem "test backtrajectory" begin
    function test_circuit(circuit)
        n=maximum((maximum(affectedqubits(g),init=1) for g in circuit),init=1)
        m=maximum((maximum(affectedbits(g),init=0) for g in circuit),init=0)

        reg1 = Register(one(Stabilizer, n), m)
        counts1 = mctrajectories(reg1, circuit; trajectories=1000, keepstates=true)

        reg2 = BacktrackRegister(n, m)
        counts2 = mctrajectories(reg2, circuit; trajectories=1000, keepstates=true)

        counts2_map = Dict()
        for ((reg2, status2), count2) in counts2
            q2 = MixedDestabilizer(quantumstate(reg2))
            counts2_map[(q2, status2)] = count2
        end

        for ((reg1, status1), count1) in counts1
            q1 = canonicalize!(quantumstate(reg1))

            key = (q1, status1)

            @assert haskey(counts2_map, key)

            count2 = counts2_map[key]

            @assert isapprox(count2, count1, rtol=0.04)
        end
    end

    @testset "circuit 1" begin
        circuit = [
            sHadamard(1),
            sCNOT(1, 2),
            sMZ(1, 1),
            sMZ(2, 2),
        ]
        test_circuit(circuit)
    end

    @testset "circuit 2" begin
        circuit = [
            sHadamard(1),
            sCNOT(1, 2),
            sCNOT(2, 3),
            sMZ(1, 1),
            sMZ(2, 2),
            sMZ(3, 3),
        ]
        test_circuit(circuit)
    end
end
