@testitem "test backtrajectory" begin
    function test_circuit(circuit)
        n=maximum((maximum(affectedqubits(g),init=1) for g in circuit),init=1)
        m=maximum((maximum(affectedbits(g),init=0) for g in circuit),init=0)

        reg1 = Register(one(Stabilizer, n), m)
        counts1 = mctrajectories(reg1, circuit; trajectories=1000, keepstates=true)

        reg2 = BacktrackRegister(n, m)
        counts2 = mctrajectories(reg2, circuit; trajectories=1000, keepstates=true)

        for ((reg1, status1), count1) in counts1
            q1 = quantumstate(reg1)
            for ((reg2, status2), count2) in counts2
                q2 = MixedDestabilizer(quantumstate(reg2))
                @assert q1 == q2
                @assert status1 == status2
                @assert isapprox(count2, count1, rtol=0.04)
            end
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

    @testset "random circuit" begin
        using Random
        Random.seed!(1234)
        n = 5
        m = 3
        len = 20
        gates = [sHadamard, sS, sCNOT, sCZ, sMX, sMY, sMZ]
        for _ in 1:3
            circuit = []
            for _ in 1:len
                gate = rand(gates)
                if gate == sCNOT || gate == sCZ
                    q1 = rand(1:n)
                    q2 = rand(1:n)
                    while q2 == q1
                        q2 = rand(1:n)
                    end
                    push!(circuit, gate(q1, q2))
                elseif gate == sHadamard || gate == sS || gate == sMX || gate == sMY || gate == sMZ
                    q = rand(1:n)
                    b = rand(0:m)
                    push!(circuit, gate(q, b))
                end
            end
            test_circuit(circuit)
        end
    end
end
