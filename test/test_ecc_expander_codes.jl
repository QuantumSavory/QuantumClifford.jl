@testitem "ECC Quantum Expander codes" begin

    using SparseArrays
    using Random
    using LinearAlgebra
    using Graphs
    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using QuantumClifford.ECC: Golay, RepCode, ReedMuller

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
        for v in G.left_nodes
            d = degree(G.g, v)
            λ[d] = get(λ, d, 0.0) + 1.0 / length(G.left_nodes)
        end
        ρ = Dict{Int,Float64}()
        for c in G.right_nodes
            d = degree(G.g, c)
            ρ[d] = get(ρ, d, 0.0) + 1.0 / length(G.right_nodes)
        end
        return λ, ρ
    end

    function generate_parity_checks(code_type::Symbol, args...)
        if code_type == :RepCode
            n = args[1]
            return sparse(parity_checks(RepCode(n)))
        elseif code_type == :ReedMuller
            r, m = args
            return parity_checks(ReedMuller(r, m))
        elseif code_type == :Golay
            n = args[1]
            if n == 23
                return Matrix{Bool}(parity_checks(Golay(n)))
            elseif n == 24
                return Matrix{Bool}(parity_checks(Golay(n)))
            else
                error("Golay code only supports n = 23 or 24")
            end
        else
            error("Unsupported code type")
        end
    end

    @testset "4.1: Validity of the Construction of Q(G1 × G2)" begin
        for n in 3:100
            H = generate_parity_checks(:RepCode, n)
            c = QuantumTannerGraphProduct(H, H)
            @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
            hx, hz = QuantumClifford.ECC.parity_checks_xz(c)
            @test QuantumClifford.ECC.verify_orthogonality(sparse(hx), sparse(hz))
        end
        for n in 3:40
            c = CyclicQuantumTannerGraphProduct(n)
            @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
            hx, hz = QuantumClifford.ECC.parity_checks_xz(c)
            @test QuantumClifford.ECC.verify_orthogonality(sparse(hx), sparse(hz))
        end
        for n in [23, 24]
            H = generate_parity_checks(:Golay, n)
            c = QuantumTannerGraphProduct(sparse(H), sparse(H))
            @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
            hx, hz = QuantumClifford.ECC.parity_checks_xz(c)
            @test QuantumClifford.ECC.verify_orthogonality(sparse(hx), sparse(hz))
        end
        for m in 3:5
            for r in 1:m-1
                H = generate_parity_checks(:ReedMuller, r, m)
                c = QuantumTannerGraphProduct(sparse(H), sparse(H))
                @test stab_looks_good(parity_checks(c); remove_redundant_rows=true)
                hx, hz = QuantumClifford.ECC.parity_checks_xz(c)
                @test QuantumClifford.ECC.verify_orthogonality(sparse(hx), sparse(hz))
            end
        end
    end

    @testset "4.2: Degree Structure of the Tanner Graphs" begin
        for n in 3:200
            H = sparse(parity_checks(RepCode(n)))
            G1 = QuantumClifford.ECC.tanner_graph_from_parity_matrix(H)
            G2 = G1 # For expander code, H1 == H2, hence G1 == G2

            # product Tanner graphs for X-type and Z-type checks
            PG_X = QuantumClifford.ECC.product_tanner_graph_X(G1, G2)
            PG_Z = QuantumClifford.ECC.product_tanner_graph_Z(G1, G2)

            # Compute degree distributions for variable and check nodes in PG_X
            var_nodes_X = PG_X.var_nodes
            var_degrees_X = compute_degree_distribution(PG_X.g, 1:length(var_nodes_X))
            check_nodes_X = PG_X.check_nodes
            check_degrees_X = compute_degree_distribution(PG_X.g, (length(var_nodes_X)+1):(length(var_nodes_X)+length(check_nodes_X)))

            # Compute degree distributions for variable and check nodes in PG_Z
            var_nodes_Z = PG_Z.var_nodes
            var_degrees_Z = compute_degree_distribution(PG_Z.g, 1:length(var_nodes_Z))
            check_nodes_Z = PG_Z.check_nodes
            check_degrees_Z = compute_degree_distribution(PG_Z.g, (length(var_nodes_Z)+1):(length(var_nodes_Z)+length(check_nodes_Z)))

            # Extract degree distributions of G1 and G2
            λ1, ρ1 = extract_degree_distribution(G1)
            λ2, ρ2 = extract_degree_distribution(G2)

            # Compute average degrees
            # λ1_avg = ∑_j λ1(j) * j
            λ1_avg = sum([k * v for (k, v) in λ1])
            # ρ1_avg = ∑_j ρ1(j) * j
            ρ1_avg = sum([k * v for (k, v) in ρ1])
            # λ2_avg = ∑_j λ2(j) * j
            λ2_avg = sum([k * v for (k, v) in λ2])
            # ρ2_avg = ∑_j ρ2(j) * j
            ρ2_avg = sum([k * v for (k, v) in ρ2])

            # Scaling factor for λ_X and λ_Z
            # scaling = (λ1_avg * λ2_avg) / (ρ1_avg * ρ2_avg)
            scaling = (λ1_avg * λ2_avg) / (ρ1_avg * ρ2_avg)

            expected_var_X = compute_expected_var_X(var_nodes_X, λ1, ρ2, scaling)
            expected_var_Z = compute_expected_var_Z(var_nodes_Z, λ2, ρ1, scaling)
            expected_check_X = compute_expected_check_X(check_nodes_X, ρ1, λ2)
            expected_check_Z = compute_expected_check_Z(check_nodes_Z, λ1, ρ2)

            # Verify distributions match
            @test var_degrees_X == expected_var_X
            @test check_degrees_X == expected_check_X
            @test var_degrees_Z == expected_var_Z
            @test check_degrees_Z == expected_check_Z
        end
    end

    @testset "Tests for parity matrix from tanner graph" begin
        Random.seed!(123)

        @testset "Empty graph" begin
            bg = QuantumClifford.ECC.generate_random_bipartite_graph(3, 2, 0.0)
            H = QuantumClifford.ECC.parity_matrix_from_tanner_graph(bg)
            @test Graphs.is_bipartite(bg.g)
            @test all(H .== false)
        end

        @testset "Fully connected bipartite graph" begin
            bg = QuantumClifford.ECC.generate_random_bipartite_graph(3, 2, 1.0)
            H = QuantumClifford.ECC.parity_matrix_from_tanner_graph(bg)
            @test Graphs.is_bipartite(bg.g)
            @test all(H .== true)
        end

        @testset "Random bipartite graphs" begin
            for _ in 1:1000
                n_v = rand(1:10)
                n_c = rand(1:10)
                p = rand() * 0.5  # Edge probability
                bg = QuantumClifford.ECC.generate_random_bipartite_graph(n_v, n_c, p)
                @test Graphs.is_bipartite(bg.g)
                H = QuantumClifford.ECC.parity_matrix_from_tanner_graph(bg)
                @test size(H) == (n_c, n_v)
            end
        end

        @testset "Large bipartite graph" begin
            bg = QuantumClifford.ECC.generate_random_bipartite_graph(1000, 50, 0.1)
            H = QuantumClifford.ECC.parity_matrix_from_tanner_graph(bg)
            @test Graphs.is_bipartite(bg.g)
            @test size(H) == (50, 1000)
            @test nnz(H) > 0
        end
    end
end
