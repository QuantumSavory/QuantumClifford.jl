using Graphs

function test_graphs()
    @testset "Graph states" begin
        for n in [1,test_sizes...]
            s = random_stabilizer(n)
            canonicalize_gott!(s)
            g, h_idx, ip_idx, z_idx = graphstate(s)
            gates = graph_gatesequence(h_idx, ip_idx, z_idx)
            gate = graph_gate(h_idx, ip_idx, z_idx, nqubits(s))
            c0 = one(CliffordOperator,nqubits(s))
            for gate in vcat(gates...) apply!(tab(c0), gate) end
            @test c0==gate
            s1 = copy(s)
            for gate in vcat(gates...) apply!(s1, gate) end
            @test canonicalize!(apply!(copy(s),c0)) == canonicalize!(s1) == canonicalize!(Stabilizer(g))
        end
    end
end

test_graphs()