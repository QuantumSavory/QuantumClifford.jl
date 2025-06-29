@testitem "Graph states" begin
    using Graphs
    using Random

    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    for n in [1,test_sizes...]
        s = random_stabilizer(n)
        g = GraphState(s)
        gates = graph_gate_sequence(g)
        gate = graph_gate(g)

        c0 = one(CliffordOperator,nqubits(s))
        for gate in gates apply!(Stabilizer(tab(c0)), gate) end
        # Test aggregated gate is indeed the same as composing individual gates from gate sequences
        @test c0==gate

        s1 = copy(s)
        for gate in gates apply!(s1, gate) end
        @test canonicalize!(apply!(s,c0)) == canonicalize!(s1) == canonicalize!(Stabilizer(g.graph))
    end

    # one more check, done manually, to ensure we are testing for cases
    # where canonicalize_gott! is causing index permutation
    # (a few of the random cases already cover that)
    s = S"- _XZ
          + _ZX
          + ZXZ"
    g = GraphState(s)
    gates = graph_gate_sequence(g)
    gate = graph_gate(g)

    c0 = one(CliffordOperator,nqubits(s))
    for gate in gates apply!(Stabilizer(tab(c0)), gate) end
    # Test aggregated gate is indeed the same as composing individual gates from gate sequences
    @test c0==gate

    s1 = copy(s)
    for gate in gates apply!(s1, gate) end
    @test canonicalize!(apply!(s,c0)) == canonicalize!(s1) == canonicalize!(Stabilizer(g.graph))
end
