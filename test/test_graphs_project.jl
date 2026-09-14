@testitem "Graph states single-qubit Pauli measurements" begin
    using Random
    import QuantumClifford.GraphSim: _project_graph!, GraphState

    test_sizes = [2,10,63,64,65,127,128,129]
    iteration = 1000

    @testset "X measurement" begin
        for n in test_sizes
            s = random_stabilizer(n)
            for _ in 1:iteration
                g = GraphState(s)
                q = rand(1:n)
                # println("state: ", s)
                # println("Graph state: ", Stabilizer(g))
                # println("qubit: ", q)
                _, res = projectXrand!(s, q)
                # println("res: ", res == 0x02)
                _project_graph!(g, q, S"X", res == 0x02)
                @test canonicalize!(Stabilizer(g)) == canonicalize!(copy(s))
            end
        end
    end
end
