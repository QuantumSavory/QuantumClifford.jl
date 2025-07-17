@testitem "Graph states" begin
    using Graphs
    using Random
    import QuantumClifford: GraphState

    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    @testset "Graph State conversion (New API)" begin
        for n in test_sizes
            s = random_stabilizer(n)
            @test canonicalize!(s) == canonicalize!(Stabilizer(GraphState(s)))
        end
    end

    @testset "Graph State conversion (Old API)" begin
        for n in test_sizes
            s = random_stabilizer(n)
            g, h_idx, ip_idx, z_idx = graphstate(s)
            gates = graph_gatesequence(h_idx, ip_idx, z_idx)
            gate = graph_gate(h_idx, ip_idx, z_idx, nqubits(s))
            c0 = one(CliffordOperator,nqubits(s))
            for gate in vcat(gates...) apply!(Stabilizer(tab(c0)), gate) end # TODO this wrapping in a Stabilizer hack is ugly; clean it up
            @test c0==gate
            s1 = copy(s)
            for gate in vcat(gates...) apply!(s1, gate) end
            @test canonicalize!(apply!(copy(s),c0)) == canonicalize!(s1) == canonicalize!(Stabilizer(g))
        end
        # one more check, done manually, to ensure we are testing for cases
        # where canonicalize_gott! is causing index permutation
        # (a few of the random cases already cover that)
        s = S"- _XZ
            + _ZX
            + ZXZ"
        g, h_idx, ip_idx, z_idx = graphstate(s)
        gates = graph_gatesequence(h_idx, ip_idx, z_idx)
        gate = graph_gate(h_idx, ip_idx, z_idx, nqubits(s))
        c0 = one(CliffordOperator,nqubits(s))
        for gate in vcat(gates...) apply!(Stabilizer(tab(c0)), gate) end # TODO this wrapping in a Stabilizer hack is ugly; clean it up
        @test c0==gate
        s1 = copy(s)
        for gate in vcat(gates...) apply!(s1, gate) end
        @test canonicalize!(apply!(copy(s),c0)) == canonicalize!(s1) == canonicalize!(Stabilizer(g))
    end
end
