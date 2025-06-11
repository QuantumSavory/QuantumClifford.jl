@testitem "ECC Quantum Expander codes" begin

    using SparseArrays
    using Random
    using LinearAlgebra
    using Graphs
    using Graphs: is_bipartite
    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using QuantumClifford.ECC: parity_checks_xz, tanner_graph_from_parity_matrix, parity_matrix_from_tanner_graph, generate_random_bipartite_graph
    using QuantumClifford.ECC: Golay, RepCode, ReedMuller

    function verify_orthogonality(HX::SparseMatrixCSC{Bool,Int}, HZ::SparseMatrixCSC{Bool,Int})
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
        return (graph=g, left=var_nodes, right=check_nodes)
    end

    function generate_parity_checks(code_type::Symbol, args...)
        if code_type == :RepCode
            n = args[1]
            return parity_checks(RepCode(n))
        elseif code_type == :ReedMuller
            r, m = args
            return sparse(parity_checks(ReedMuller(r, m)))
        elseif code_type == :Golay
            n = args[1]
            if n == 23
                return sparse(Matrix{Bool}(parity_checks(Golay(n))))
            elseif n == 24
                return sparse(Matrix{Bool}(parity_checks(Golay(n))))
            else
                error("Golay code only supports n = 23 or 24")
            end
        else
            error("Unsupported code type")
        end
    end

    @testset "4.1: Validity of the Construction of Q(G1 × G2)" begin
        @testset "Repetition codes" begin
            for n in 3:20
                H = sparse(parity_checks(RepCode(n)))
                c = QuantumTannerGraphProduct(H, H)
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = parity_checks_xz(c)
                @test verify_orthogonality(sparse(hx), sparse(hz))
            end
        end

        @testset "Cyclic codes" begin
            for n in 3:20
                c = CyclicQuantumTannerGraphProduct(n)
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = parity_checks_xz(c)
                @test verify_orthogonality(sparse(hx), sparse(hz))
            end
        end

        @testset "Golay codes" begin
            for n in [23, 24]
                H = sparse(Matrix{Bool}(parity_checks(Golay(n))))
                c = QuantumTannerGraphProduct(H, H)
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = parity_checks_xz(c)
                @test verify_orthogonality(sparse(hx), sparse(hz))
            end
        end

        @testset "Reed-Muller codes" begin
            for m in 3:5, r in 1:m-1
                H = sparse(parity_checks(ReedMuller(r, m)))
                c = QuantumTannerGraphProduct(H, H)
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = parity_checks_xz(c)
                @test verify_orthogonality(sparse(hx), sparse(hz))
            end
        end
    end

    @testset "Fully connected bipartite graph" begin
        bg = generate_random_bipartite_graph(3, 2, 1.0)
        H = parity_matrix_from_tanner_graph(bg.graph, bg.left, bg.right)
        @test Graphs.is_bipartite(bg.graph)
        @test all(H .== true)
    end

    @testset "Random bipartite graphs" begin
        for _ in 1:1000
            n_v = rand(1:10)
            n_c = rand(1:10)
            p = rand() * 0.5
            bg = generate_random_bipartite_graph(n_v, n_c, p)
            @test is_bipartite(bg.graph)
            H = parity_matrix_from_tanner_graph(bg.graph, bg.left, bg.right)
            @test size(H) == (n_c, n_v)
        end
    end

    @testset "Large bipartite graph" begin
        bg = generate_random_bipartite_graph(1000, 50, 0.1)
        H = parity_matrix_from_tanner_graph(bg.graph, bg.left, bg.right)
        @test is_bipartite(bg.graph)
        @test size(H) == (50, 1000)
        @test nnz(H) > 0
    end
end
