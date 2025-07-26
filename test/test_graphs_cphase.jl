@testitem "Graph state CPHASE simulation" begin
    using Random
    import QuantumClifford: GraphState, local_comp!

    test_sizes = [2,10,63,64,65,127,128,129]

    @testset "CPHASE correctness" begin
        for n in test_sizes
            s = random_stabilizer(n)
            g = GraphState(s)
            a, b = rand(1:n), rand(1:n)
            while b == a
                b = rand(1:n)
            end
            gate = sCPHASE(a, b)
            @test canonicalize!(Stabilizer(apply!(g, gate))) == canonicalize!(apply!(s, gate))
        end

        # test a special case where q1, q2 are disconnected and both VOPs not in Z_COMMUTATION_SUBGROUP
        s = S"ZXII XZII IIZX IIXZ"
        gs = GraphState(s)
        gate = sCPHASE(2,3)
        apply!(gs, gate)
        apply!(s, gate)
        @test canonicalize!(Stabilizer(gs)) == canonicalize!(s)
    end

    @testset "Local Complementation equivalence" begin
        for n in test_sizes
            s = random_stabilizer(n)
            g = GraphState(s)
            @test canonicalize!(Stabilizer(local_comp!(g, rand(1:n)))) == canonicalize!(s)
        end
    end
end
