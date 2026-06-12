@testitem "Adapters" tags=[:ecc, :ecc_adapters] begin
    using QuantumClifford
    using QuantumClifford.ECC: CSS, Surface, Steane7,
                                BivariateBicycleViaCirculantMat,
                                two_block_group_algebra_code,
                                CodePair, build_adapter, build_adapter_intercode,
                                build_adapter_intracode, deform_code, skiptree,
                                AuxiliaryGraph, SkipTreeOutput, Adapter,
                                code_n, code_k, code_s, parity_matrix_x,
                                parity_matrix_z, logz_ops, logx_ops,
                                distance, DistanceMIPAlgorithm
    using QuantumClifford.ECC.Adapters: build_initial_aux_graph,
                                         cellulate_long_cycles!,
                                         relative_expansion, verify_desideratum_4,
                                         canonical_H_R, port_function_as_matrix,
                                         incidence_matrix_transpose,
                                         stabilizer_modification_matrix,
                                         adapter_bridge_x_check_matrix,
                                         edge_pair_to_index
    using QuantumClifford: stab_to_gf2, nqubits
    using Hecke: group_algebra, GF, abelian_group, gens
    using Graphs: SimpleGraph, add_edge!, edges, ne, nv, src, dst, has_edge,
                  path_graph, cycle_graph, complete_graph, is_connected,
                  erdos_renyi, cycle_basis, degree, neighbors
    using SparseArrays
    using LinearAlgebra
    using Random
    using JuMP, HiGHS

    as_css(c) = CSS(Matrix{Bool}(parity_matrix_x(c)),
                    Matrix{Bool}(parity_matrix_z(c)))

    css_commutes(m::CSS) =
        count(!iszero, (Int.(m.Hx) * transpose(Int.(m.Hz))) .% 2) == 0

    # 2BGA bivariate-bicycle code from polynomial term lists `(:x|:y, power)`.
    function bb_2bga(l, m, Aterms, Bterms)
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        var = Dict(:x => x, :y => y)
        Apoly = sum(var[s]^p for (s, p) in Aterms)
        Bpoly = sum(var[s]^p for (s, p) in Bterms)
        two_block_group_algebra_code(Apoly, Bpoly)
    end

    # Edge-vertex incidence (ne × nv) of a SimpleGraph in `edges(g)` order.
    function incidence_matrix(g::SimpleGraph)
        I = Int[]; J = Int[]
        for (k, e) in enumerate(edges(g))
            push!(I, k); push!(J, src(e))
            push!(I, k); push!(J, dst(e))
        end
        sparse(I, J, trues(length(I)), ne(g), nv(g))
    end

    is_permutation(P::SparseMatrixCSC) =
        all(sum(P; dims = 1) .== 1) && all(sum(P; dims = 2) .== 1)

    # Run the full SkipTree contract on `g` and return the output.
    function check_skiptree(g::SimpleGraph; root = 1)
        n = nv(g); m = ne(g)
        out = skiptree(g; root = root)
        @test out isa SkipTreeOutput
        @test out.target === :H_R
        @test size(out.T) == (n - 1, m)
        @test size(out.P) == (n, n)
        @test maximum(vec(sum(out.T; dims = 2))) ≤ 3   # row weight ≤ 3
        @test maximum(vec(sum(out.T; dims = 1))) ≤ 2   # col weight ≤ 2
        @test is_permutation(out.P)
        G = incidence_matrix(g)
        TGP = mod.(Matrix{Int}(out.T) * Matrix{Int}(G) * Matrix{Int}(out.P), 2)
        @test TGP == Matrix{Int}(canonical_H_R(n))
        out
    end

    @testset "SkipTree on synthetic graphs" begin
        check_skiptree(SimpleGraph(path_graph(5)))
        check_skiptree(SimpleGraph(cycle_graph(7)))
        Random.seed!(2026)
        for trial in 1:5
            g = SimpleGraph(erdos_renyi(20, 0.3; seed = 17 * trial))
            if !is_connected(g)
                for v in 2:nv(g); add_edge!(g, v - 1, v); end
            end
            check_skiptree(g)
        end
    end

    @testset "Auxiliary graph (build_initial_aux_graph)" begin
        # Surface(3,3) gives a weight-3 Z̄.
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        @test length(z) == 3
        aux = build_initial_aux_graph(z, as_css(Surface(3, 3)))
        @test nv(aux.graph) == 3
        @test aux.logical_support == z
        @test all(length(p) == 1 for p in aux.port_function)    # injective
        # Every stored matching pair is an edge of the current graph.
        for matchings in aux.stabilizer_matchings, p in matchings
            @test has_edge(aux.graph, p[1], p[2])
        end
        # Boundary case: requires ≥ 2 support qubits.
        @test_throws ArgumentError build_initial_aux_graph([1], as_css(Surface(3, 3)))
    end

    @testset "cellulate_long_cycles!" begin
        c = bb_2bga(6, 6, [(:x, 3), (:y, 1), (:y, 2)],
                          [(:y, 3), (:x, 1), (:x, 2)])
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(c))[1, code_n(c)+1:2*code_n(c)]))
        aux0 = build_initial_aux_graph(z, as_css(c))
        aux  = cellulate_long_cycles!(aux0; max_cycle_len = 6)

        # Every basis cycle is bounded by max_cycle_len.
        @test all(length(cyc) ≤ 6 for cyc in cycle_basis(aux.graph))

        # Algebraic invariant: N · G ≡ 0 (mod 2), where G is the
        # ne × nv incidence (each cycle is a closed walk).
        G = incidence_matrix(aux.graph)
        N = aux.cycle_basis_matrix
        prod = mod.(Matrix{Int}(N) * Matrix{Int}(G), 2)
        @test all(iszero, prod)

        # Explicit chord routes through to the graph.
        aux_chord = cellulate_long_cycles!(build_initial_aux_graph(z, as_css(c));
                                            max_cycle_len = 6,
                                            chords = [(1, 3)])
        @test has_edge(aux_chord.graph, 1, 3)
    end

    @testset "relative_expansion (Definition 2)" begin
        # Path P_4 against the whole vertex set: β_4 = 1/2.
        @test relative_expansion(SimpleGraph(path_graph(4)), [1, 2, 3, 4], 4) ==
              1 // 2
        # Complete K_5 against the whole vertex set: β_5 = 3.
        @test relative_expansion(SimpleGraph(complete_graph(5)),
                                  [1, 2, 3, 4, 5], 5) == 3 // 1
        # Size guard: brute force only up to n ≤ 24.
        @test_throws ArgumentError relative_expansion(SimpleGraph(path_graph(25)),
                                                       [1, 2], 1)
    end

    @testset "Merge helpers" begin
        # canonical_H_R(n) has the banded `1 1` pattern on `n-1` rows.
        for n in (2, 3, 5, 8)
            H = canonical_H_R(n)
            @test size(H) == (n - 1, n)
            for i in 1:(n - 1)
                @test H[i, i] && H[i, i + 1]
                @test sum(H[i, :]) == 2
            end
        end

        # adapter_bridge_x_check_matrix on a synthetic SkipTree pair.
        g = SimpleGraph(cycle_graph(4))
        st = skiptree(g)
        bridge = adapter_bridge_x_check_matrix(st, st, 4)
        @test size(bridge) == (3, size(st.T, 2) + 4 + size(st.T, 2))
        @test_throws ArgumentError adapter_bridge_x_check_matrix(st, st, 1)
    end

    @testset "edge_pair_to_index" begin
        g = SimpleGraph(cycle_graph(5))
        # Iteration order matches what edge_pair_to_index returns.
        for (k, e) in enumerate(edges(g))
            u, v = minmax(src(e), dst(e))
            @test edge_pair_to_index(g, (u, v)) == k
        end
        # Non-edge → 0.
        @test edge_pair_to_index(SimpleGraph(path_graph(4)), (1, 4)) == 0
        # Argument check: pair must be sorted.
        @test_throws ArgumentError edge_pair_to_index(g, (3, 1))
    end

    @testset "Surface(3,3) × Surface(3,3) inter-code → [[33, 1, 3]]" begin
        c1 = as_css(Surface(3, 3))
        c2 = as_css(Surface(3, 3))
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        adapter = build_adapter(CodePair(c1, c2, z, z))
        @test code_n(adapter) == 33                  # 13 + 2 + 3 + 2 + 13
        @test code_k(adapter) == 1                   # one logical absorbed
        @test css_commutes(adapter.merged_code)
        @test adapter.adapter_width == 3
        # Distance preserved (Surface(3,3) has d = 3).
        @test distance(adapter, DistanceMIPAlgorithm(solver = HiGHS)) == 3
    end

    @testset "BB[[72, 12, 6]] intra-code joint Z measurement" begin
        c_orig = bb_2bga(6, 6, [(:x, 3), (:y, 1), (:y, 2)],
                              [(:y, 3), (:x, 1), (:x, 2)])
        @test (code_n(c_orig), code_k(c_orig)) == (72, 12)
        c = as_css(c_orig)
        lz = stab_to_gf2(logz_ops(c_orig))
        n0 = code_n(c_orig)
        # Pick two of the BB's twelve Z-logicals to jointly measure.
        z1 = sort(findall(!iszero, lz[1, n0 + 1:2n0]))
        z2 = sort(findall(!iszero, lz[2, n0 + 1:2n0]))
        adapter = build_adapter(CodePair(c, c, z1, z2))
        @test css_commutes(adapter.merged_code)
        @test code_k(adapter) == code_k(c_orig) - 1    # joint measurement consumes 1
        # Weight bounds from paper Theorem 5 / Table III on BB family.
        @test maximum(sum(parity_matrix_x(adapter), dims = 2)) ≤ 8
        @test maximum(sum(parity_matrix_z(adapter), dims = 2)) ≤ 8
    end

    @testset "Gross [[144, 12, 12]] BB structural regression" begin
        # Bravyi-et-al-2024 gross BB; library-constructed via the circulant
        # matrix path. We do not run the full distance MIP here (hours at d=12);
        # only check the construction is well-formed and the joint measurement
        # bookkeeping is correct.
        l, m_param = 12, 6
        A = [(:x, 3), (:y, 1), (:y, 2)]
        B = [(:y, 3), (:x, 1), (:x, 2)]
        g_orig = BivariateBicycleViaCirculantMat(l, m_param, A, B)
        @test (code_n(g_orig), code_k(g_orig)) == (144, 12)

        g = as_css(g_orig)
        lz = stab_to_gf2(logz_ops(g_orig))
        n0 = code_n(g_orig)
        z1 = sort(findall(!iszero, lz[5,  n0 + 1:2n0]))
        z2 = sort(findall(!iszero, lz[10, n0 + 1:2n0]))
        adapter = build_adapter(CodePair(g, g, z1, z2))

        # Column count from intra-code block layout.
        expected_n = n0 + ne(adapter.aux_l.graph) +
                          adapter.adapter_width +
                          ne(adapter.aux_r.graph)
        @test code_n(adapter) == expected_n
        @test adapter.adapter_width == min(length(z1), length(z2))
        @test css_commutes(adapter.merged_code)
        @test code_k(adapter) == code_k(g_orig) - 1
        # Weight bounds.
        Hx = parity_matrix_x(adapter)
        Hz = parity_matrix_z(adapter)
        @test maximum(sum(Hx, dims = 2)) ≤ 8
        @test maximum(sum(Hz, dims = 2)) ≤ 8
        @test maximum(sum(vcat(Hx, Hz), dims = 1)) ≤ 9
    end

    @testset "deform_code (single-logical) on Surface(3,3)" begin
        c = as_css(Surface(3, 3))
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        dc = deform_code(c, z)
        # n_deformed = n_code + ne(aux.graph) for the cellulated aux graph.
        aux = cellulate_long_cycles!(build_initial_aux_graph(z, c);
                                      max_cycle_len = Int(maximum(sum(c.Hx; dims = 2))))
        @test code_n(dc) == code_n(c) + ne(aux.graph)
        @test code_k(dc) == 0    # measuring the only logical consumes it
        @test css_commutes(dc)
    end

    @testset "build_adapter dispatcher" begin
        c  = as_css(Surface(3, 3))
        c1 = as_css(Surface(3, 3))
        c2 = as_css(Surface(3, 3))
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        # Distinct objects route to the inter-code path.
        @test code_n(build_adapter(CodePair(c1, c2, z, z))) == 33
        # Same object routes to the intra-code path.
        @test build_adapter(CodePair(c, c, z, z)).adapter_width == 3
        # Inter-code path rejects same-object pairs.
        @test_throws ArgumentError build_adapter_intercode(CodePair(c, c, z, z))
        # Intra-code path rejects distinct-object pairs.
        @test_throws ArgumentError build_adapter_intracode(CodePair(c1, c2, z, z))
    end

    @testset "Adapter plugs directly into AbstractCSSCode dispatch" begin
        c1 = as_css(Surface(3, 3))
        c2 = as_css(Surface(3, 3))
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        adapter = build_adapter(CodePair(c1, c2, z, z))
        @test adapter isa QuantumClifford.ECC.AbstractCSSCode
        @test code_n(adapter) == 33
        @test code_k(adapter) == 1
        @test code_s(adapter) == 32
        @test parity_matrix_x(adapter) === adapter.merged_code.Hx
        @test parity_matrix_z(adapter) === adapter.merged_code.Hz
    end

    @testset "Distance MIP plumbing (sanity)" begin
        # If the MIP formulation regresses, these will visibly break.
        @test distance(Steane7(), DistanceMIPAlgorithm(solver = HiGHS)) == 3
        @test distance(Surface(3, 3), DistanceMIPAlgorithm(solver = HiGHS)) == 3
    end
end
