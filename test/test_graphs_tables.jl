@testitem "Graph states precompiled tables correctness" begin
    import QuantumClifford: GraphState, IP_SQRTX_DECOMPOSITION_TABLE, ISOLATED_CPHASE_TABLE, Z_COMMUTATION_SUBGROUP, SINGLE_QUBIT_MULTIPLICATION_TABLE, gen_graph_state, vops

    @testset "Single Clifford operator decomposition table" begin
        for (k, v) in IP_SQRTX_DECOMPOSITION_TABLE
            # check if each sequence "adds up" to the operator
            t = one(CliffordOperator, 1)
            for u in v
                t = u * t
            end
            @test SingleQubitOperator(t) == k
        end
    end

    @testset "Single qubit operator multiplication table" begin
        for (k, v) in SINGLE_QUBIT_MULTIPLICATION_TABLE
            (v1, v2) = k
            @test v1 * CliffordOperator(v2, 1) == CliffordOperator(v, 1)
        end
    end

    @testset "Isolated CPHASE Table" begin
        for ((connected, U1, U2), res) in ISOLATED_CPHASE_TABLE
            g_init = gen_graph_state(connected, U1, U2)
            # check every entry satisfies the constraint
            for qubit_idx in 1:2
                if CliffordOperator(vops(g_init)[qubit_idx], 1) in Z_COMMUTATION_SUBGROUP
                    @test CliffordOperator(vops(res)[qubit_idx], 1) in Z_COMMUTATION_SUBGROUP
                end
            end
            # check every entry agrees with stabilizer result
            @test canonicalize!(apply!(Stabilizer(g_init), sCPHASE(1, 2))) == canonicalize!(Stabilizer(res))
        end
    end
end
