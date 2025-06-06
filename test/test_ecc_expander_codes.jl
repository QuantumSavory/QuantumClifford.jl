@testitem "ECC Quantum Expander codes" begin

    using SparseArrays
    using Random
    using LinearAlgebra
    using Graphs
    using Graphs: is_bipartite
    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using QuantumClifford.ECC: parity_checks_xz, product_tanner_graph_X, product_tanner_graph_Z, tanner_graph_from_parity_matrix, parity_matrix_from_product_tanner_graph, parity_matrix_from_tanner_graph, generate_random_bipartite_graph
    using QuantumClifford.ECC: Golay, RepCode, ReedMuller

    # Checks whether the parity-check matrices H_X and H_Z are orthogonal in GF(2),
    # a key requirement for CSS code construction. In product Tanner graphs, the
    # existence of 4-cycles ensures orthogonality.
    function verify_orthogonality(HX::SparseMatrixCSC{Bool,Int}, HZ::SparseMatrixCSC{Bool,Int})
        product_mod2 = mod.(Int.(HX) * transpose(Int.(HZ)), 2)
        return all(product_mod2 .== 0)
    end

    # Function to compute the degree distribution of nodes in a graph
    function compute_degree_distribution(graph, nodes)
        degrees = Dict{Int, Int}()
        for node in nodes
            deg = degree(graph, node)
            degrees[deg] = get(degrees, deg, 0) + 1
        end
        return degrees
    end

    # Compute the expected variable degrees for X-type (λ_X).
    # λ_X(j) = λ1(j) + (λ1_avg * λ2_avg / ρ1_avg * ρ2_avg) * ρ2(j) / (λ1_avg * λ2_avg / ρ1_avg * ρ2_avg + 1)
    function compute_expected_var_X(var_nodes_X, λ1, ρ2, scaling)
        expected_var_X = Dict{Int,Int}()
        total_vars_X = length(var_nodes_X)
        for j in union(keys(λ1), keys(ρ2))
            numerator = get(λ1, j, 0.0) + scaling * get(ρ2, j, 0.0)
            fraction = numerator / (scaling + 1.0)
            expected_count = round(Int, total_vars_X * fraction)
            if expected_count > 0
                expected_var_X[j] = expected_count
            end
        end
        return expected_var_X
    end

    # Compute the expected check degrees for X-type (ρ_X).
    # ρ_X(k) = ∑_{i,j:i+j=k} ρ1(i) * λ2(j)
    function compute_expected_check_X(check_nodes_X, ρ1, λ2)
        expected_check_X = Dict{Int,Int}()
        total_checks_X = length(check_nodes_X)
        for k in keys(ρ1), l in keys(λ2)
            if !haskey(expected_check_X, k + l)
                expected_check_X[k + l] = 0
            end
            expected_check_X[k + l] += round(Int, total_checks_X * ρ1[k] * λ2[l])
        end
        return expected_check_X
    end

    # Compute the expected variable degrees for Z-type (λ_Z).
    # λ_Z(j) = λ2(j) + (λ1_avg * λ2_avg / ρ1_avg * ρ2_avg) * ρ1(j) / (λ1_avg * λ2_avg / ρ1_avg * ρ2_avg + 1)
    function compute_expected_var_Z(var_nodes_Z, λ2, ρ1, scaling)
        expected_var_Z = Dict{Int,Int}()
        total_vars_Z = length(var_nodes_Z)
        for j in union(keys(λ2), keys(ρ1))
            numerator = get(λ2, j, 0.0) + scaling * get(ρ1, j, 0.0)
            fraction = numerator / (scaling + 1.0)
            expected_count = round(Int, total_vars_Z * fraction)
            if expected_count > 0
                expected_var_Z[j] = expected_count
            end
        end
        return expected_var_Z
    end

    # Compute the expected check degrees for Z-type (ρ_Z).
    # ρ_Z(k) = ∑_{i,j:i+j=k} λ1(i) * ρ2(j)
    function compute_expected_check_Z(check_nodes_Z, λ1, ρ2)
        expected_check_Z = Dict{Int,Int}()
        total_checks_Z = length(check_nodes_Z)
        for k in keys(λ1), l in keys(ρ2)
            if !haskey(expected_check_Z, k + l)
                expected_check_Z[k + l] = 0
            end
            expected_check_Z[k + l] += round(Int, total_checks_Z * λ1[k] * ρ2[l])
        end
        return expected_check_Z
    end

    # Extract the degree distributions for variable nodes (λ) and check nodes (ρ).
    function extract_degree_distribution(G)
        λ = Dict{Int,Float64}()
        for v in G.left
            d = degree(G.graph, v)
            λ[d] = get(λ, d, 0.0) + 1.0 / length(G.left)
        end
        ρ = Dict{Int,Float64}()
        for c in G.right
            d = degree(G.graph, c)
            ρ[d] = get(ρ, d, 0.0) + 1.0 / length(G.right)
        end
        return λ, ρ
    end

    function compare_degree_dicts(actual::Dict{Int,Int}, expected::Dict{Int,Int}, atol::Int=1)
        all_degrees = union(keys(actual), keys(expected))
        for d in all_degrees
            a = get(actual, d, 0)
            e = get(expected, d, 0)
            if !(abs(a - e) ≤ atol)
                return false
            end
        end
        return true
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

    @testset "4.2: Degree Structure of the Tanner Graphs" begin
        Random.seed!(123)
        @testset "Repetition codes" begin
            for n in [3, 5, 10, 20]
                H = sparse(parity_checks(RepCode(n)))
                G1 = tanner_graph_from_parity_matrix(H)
                G2 = G1

                PG_X = product_tanner_graph_X(G1, G2)
                PG_Z = product_tanner_graph_Z(G1, G2)

                var_degrees_X = compute_degree_distribution(PG_X.graph, 1:length(PG_X.var_nodes))
                check_degrees_X = compute_degree_distribution(PG_X.graph, (length(PG_X.var_nodes)+1):nv(PG_X.graph))

                var_degrees_Z = compute_degree_distribution(PG_Z.graph, 1:length(PG_Z.var_nodes))
                check_degrees_Z = compute_degree_distribution(PG_Z.graph, (length(PG_Z.var_nodes)+1):nv(PG_Z.graph))

                λ1, ρ1 = extract_degree_distribution(G1)
                λ2, ρ2 = extract_degree_distribution(G2)
                
                λ1_avg = sum(k * v for (k, v) in λ1)
                ρ1_avg = sum(k * v for (k, v) in ρ1)
                scaling = (λ1_avg * λ1_avg) / (ρ1_avg * ρ1_avg)  # G1 == G2

                expected_var_X = compute_expected_var_X(PG_X.var_nodes, λ1, ρ1, scaling)
                expected_var_Z = compute_expected_var_Z(PG_Z.var_nodes, λ1, ρ1, scaling)
                expected_check_X = compute_expected_check_X(PG_X.check_nodes, ρ1, λ1)
                expected_check_Z = compute_expected_check_Z(PG_Z.check_nodes, λ1, ρ1)

                @test compare_degree_dicts(var_degrees_X, expected_var_X, 1)
                @test compare_degree_dicts(check_degrees_X, expected_check_X, 1)
                @test compare_degree_dicts(var_degrees_Z, expected_var_Z, 1)
                @test compare_degree_dicts(check_degrees_Z, expected_check_Z, 1)
            end
        end
    end

    @testset "Parity matrix conversions" begin
        @testset "Product graph conversions" begin
            H1 = sparse(parity_checks(RepCode(3)))
            H2 = sparse(parity_checks(RepCode(4)))
            G1 = tanner_graph_from_parity_matrix(H1)
            G2 = tanner_graph_from_parity_matrix(H2)
            
            PG_X = product_tanner_graph_X(G1, G2)
            HX = parity_matrix_from_product_tanner_graph(PG_X.graph, PG_X.var_nodes, PG_X.check_nodes)
            @test size(HX) == (length(PG_X.check_nodes), length(PG_X.var_nodes))
            
            PG_Z = product_tanner_graph_Z(G1, G2)
            HZ = parity_matrix_from_product_tanner_graph(PG_Z.graph, PG_Z.var_nodes, PG_Z.check_nodes)
            @test size(HZ) == (length(PG_Z.check_nodes), length(PG_Z.var_nodes))
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
