@testitem "Quantum Tanner Graph Codes" begin

    using Random
    using Graphs
    using Graphs: is_bipartite
    using SparseArrays
    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using QuantumClifford.ECC: parity_matrix_x, parity_matrix_z
    using QECCore
    using QECCore: tanner_graph_from_parity_matrix, parity_matrix_from_tanner_graph, BipartiteGraph

    function verify_orthogonality(HX::AbstractMatrix, HZ::AbstractMatrix)
        product_mod2 = mod.(Int.(HX) * transpose(Int.(HZ)), 2)
        return all(product_mod2 .== 0)
    end

    function generate_random_bipartite_graph(n_vars::Int, n_checks::Int, edge_prob::Float64)
        g = SimpleGraph(n_vars + n_checks)
        var_nodes = collect(1:n_vars)
        check_nodes = collect(n_vars+1:n_vars+n_checks)
        for v in var_nodes, c in check_nodes
            rand() < edge_prob && add_edge!(g, v, c)
        end
        return BipartiteGraph(g, var_nodes, check_nodes)
    end

    function generate_parity_checks(code_type::Symbol, args...)
        if code_type == :RepCode
            n = args[1]
            return parity_matrix(RepCode(n))
        elseif code_type == :ReedMuller
            r, m = args
            return ReedMuller(r, m)
        elseif code_type == :Golay
            n = args[1]
            if n == 23
                return parity_matrix(Golay(n))
            elseif n == 24
                return parity_matrix(Golay(n))
            else
                error("Golay code only supports n = 23 or 24")
            end
        else
            error("Unsupported code type")
        end
    end

    @testset "4.1: Validity of the Construction of Q(G1 × G2)" begin
        @testset "Repetition codes" begin
            for n in 3:10
                H = parity_matrix(RepCode(n))
                c = QuantumTannerGraphProduct(H, H)
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = parity_matrix_x(c), parity_matrix_z(c)
                @test verify_orthogonality(hx, hz)
            end
        end

        @testset "Cyclic codes" begin
            for n in 3:10
                c = CyclicQuantumTannerGraphProduct(n)
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = parity_matrix_x(c), parity_matrix_z(c)
                @test verify_orthogonality(hx, hz)
            end
        end

        @testset "Golay codes" begin
            for n in [23, 24]
                H = parity_matrix(Golay(n))
                c = QuantumTannerGraphProduct(H, H)
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = parity_matrix_x(c), parity_matrix_z(c)
                @test verify_orthogonality(hx, hz)
            end
        end

        @testset "Reed-Muller codes" begin
            for m in 3:5, r in 1:m-1
                H = parity_matrix(ReedMuller(r, m))
                c = QuantumTannerGraphProduct(H, H)
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = parity_matrix_x(c), parity_matrix_z(c)
                @test verify_orthogonality(hx, hz)
            end
        end
    end

    @testset "Fully connected bipartite graph" begin
        bg = generate_random_bipartite_graph(3, 2, 1.0)
        H = parity_matrix_from_tanner_graph(bg)
        @test Graphs.is_bipartite(bg.g)
        @test all(H .== true)
    end

    @testset "Random bipartite graphs" begin
        for _ in 1:1000
            n_v = rand(1:10)
            n_c = rand(1:10)
            p = rand() * 0.5
            bg = generate_random_bipartite_graph(n_v, n_c, p)
            @test is_bipartite(bg.g)
            H = parity_matrix_from_tanner_graph(bg)
            @test size(H) == (n_c, n_v)
        end
    end

    @testset "Large bipartite graph" begin
        bg = generate_random_bipartite_graph(1000, 50, 0.1)
        H = parity_matrix_from_tanner_graph(bg)
        @test is_bipartite(bg.g)
        @test size(H) == (50, 1000)
        @test nnz(H) > 0
    end
end
