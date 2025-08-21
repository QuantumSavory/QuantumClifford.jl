@testitem "Graph state CPHASE simulation" begin
    using Random
    import QuantumClifford: GraphState, local_comp!

    test_sizes = [2,10,63,64,65,127,128,129]
    iteration = 10000

    @testset "CPHASE correctness" begin
        for _ in 1:iteration
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
        end

        # test a special case where q1, q2 are disconnected and both VOPs not in Z_COMMUTATION_SUBGROUP
        s = S"ZXII XZII IIZX IIXZ"
        gs = GraphState(s)
        gate = sCPHASE(2,3)
        apply!(gs, gate)
        apply!(s, gate)
        @test canonicalize!(Stabilizer(gs)) == canonicalize!(s)

        # test when q1, q2 are disconnected but each connected to other nodes as well
        s = S"+ __XXXYZZ__
              + YY_ZZY_Z_X
              + XZZZY_YZZY
              - YYYY_ZZXXZ
              + XZYZXXZXXZ
              + ZZXYXYX_Y_
              + ___YXZYXY_
              + __ZXYZZXX_
              + Y_YY_____Z
              + ZZY__XXXYX"
        gs = GraphState(s)
        gate = sCPHASE(6,8)
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
