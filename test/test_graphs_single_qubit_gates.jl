@testitem "Graph states single qubit gates" begin
    using Random
    import QuantumClifford: GraphState

    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    @testset "Single qubit gate application" begin
        for n in test_sizes
            for i in 1:6
                op = enumerate_single_qubit_gates(i, qubit=rand(1:n), phases=(rand(Bool),rand(Bool)))
                s = random_stabilizer(n)
                gs = GraphState(s)
                @test canonicalize!(apply!(copy(s),op))==canonicalize!(Stabilizer(apply!(copy(gs),op)))
            end
        end
    end
end
