@testitem "Graph states" begin
    using Graphs
    using Random
    import QuantumClifford.GraphSim: GraphState

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
            @test canonicalize!(apply!(Stabilizer(g), inv(gate))) == canonicalize!(copy(s))
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

    @testset "Graph state reconstruction direction and order" begin
        # Y needs Phase on reconstruction. For -Z, Pauli Z must precede Hadamard.
        states = (S"Y", S"-Y", S"-Z", S" Y_
                                      -_Z",
                  # The [[4,2,2]] code with logical Y and X fixed to +1.
                  S"XXXX
                    ZZZZ
                    ZXY_
                    _X_X")
        for s in states
            g, h_idx, ip_idx, z_idx = graphstate(s)
            to_graph = graph_gate(h_idx, ip_idx, z_idx, nqubits(s))
            @test canonicalize!(apply!(copy(s), to_graph)) == canonicalize!(Stabilizer(g))
            @test canonicalize!(apply!(Stabilizer(g), inv(to_graph))) == canonicalize!(copy(s))

            recovered = Stabilizer(g)
            gates = vcat(graph_gatesequence(h_idx, ip_idx, z_idx)...)
            for gate in reverse(gates)
                apply!(recovered, inv(gate))
            end
            @test canonicalize!(recovered) == canonicalize!(copy(s))
        end
    end
end
