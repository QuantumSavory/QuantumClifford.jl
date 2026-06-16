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
                                         assemble_merged_code_intercode,
                                         assemble_merged_code_intracode,
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

    function bb_2bga(l, m, Aterms, Bterms)
        GA = group_algebra(GF(2), abelian_group([l, m]))
        x, y = gens(GA)
        var = Dict(:x => x, :y => y)
        Apoly = sum(var[s]^p for (s, p) in Aterms)
        Bpoly = sum(var[s]^p for (s, p) in Bterms)
        two_block_group_algebra_code(Apoly, Bpoly)
    end

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

    function check_skiptree(g::SimpleGraph; root = 1)
        n = nv(g); m = ne(g)
        out = skiptree(g; root = root)
        @test out isa SkipTreeOutput
        @test out.target === :H_R
        @test size(out.T) == (n - 1, m)
        @test size(out.P) == (n, n)
        @test maximum(vec(sum(out.T; dims = 2))) ≤ 3
        @test maximum(vec(sum(out.T; dims = 1))) ≤ 2
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
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        @test length(z) == 3
        aux = build_initial_aux_graph(z, as_css(Surface(3, 3)))
        @test nv(aux.graph) == 3
        @test aux.logical_support == z
        @test all(length(p) == 1 for p in aux.port_function)
        for matchings in aux.stabilizer_matchings, p in matchings
            @test has_edge(aux.graph, p[1], p[2])
        end
        @test_throws ArgumentError build_initial_aux_graph([1], as_css(Surface(3, 3)))
    end

    @testset "cellulate_long_cycles!" begin
        c = bb_2bga(6, 6, [(:x, 3), (:y, 1), (:y, 2)],
                          [(:y, 3), (:x, 1), (:x, 2)])
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(c))[1, code_n(c)+1:2*code_n(c)]))
        aux0 = build_initial_aux_graph(z, as_css(c))
        aux  = cellulate_long_cycles!(aux0; max_cycle_len = 6)

        @test all(length(cyc) ≤ 6 for cyc in cycle_basis(aux.graph))
        # N · G ≡ 0 (mod 2)
        G = incidence_matrix(aux.graph)
        N = aux.cycle_basis_matrix
        @test all(iszero, mod.(Matrix{Int}(N) * Matrix{Int}(G), 2))

        aux_chord = cellulate_long_cycles!(build_initial_aux_graph(z, as_css(c));
                                            max_cycle_len = 6,
                                            chords = [(1, 3)])
        @test has_edge(aux_chord.graph, 1, 3)

        @test_throws ArgumentError cellulate_long_cycles!(
            build_initial_aux_graph(z, as_css(c)); max_cycle_len = 6,
            chords = [(0, 2)])
        @test_throws ArgumentError cellulate_long_cycles!(
            build_initial_aux_graph(z, as_css(c)); max_cycle_len = 6,
            chords = [(1, 10_000)])
        @test_throws ArgumentError cellulate_long_cycles!(
            build_initial_aux_graph(z, as_css(c)); max_cycle_len = 6,
            chords = [(2, 2)])
        aux_existing = build_initial_aux_graph(z, as_css(c))
        e_existing   = first(edges(aux_existing.graph))
        u_e, v_e     = minmax(src(e_existing), dst(e_existing))
        @test_throws ArgumentError cellulate_long_cycles!(
            build_initial_aux_graph(z, as_css(c)); max_cycle_len = 6,
            chords = [(u_e, v_e)])
        @test_throws ArgumentError cellulate_long_cycles!(
            build_initial_aux_graph(z, as_css(c)); max_cycle_len = 2)
    end

    @testset "relative_expansion (Definition 2)" begin
        # β_4(P_4) = 1/2
        @test relative_expansion(SimpleGraph(path_graph(4)), [1, 2, 3, 4], 4) ==
              1 // 2
        # β_5(K_5) = 3
        @test relative_expansion(SimpleGraph(complete_graph(5)),
                                  [1, 2, 3, 4, 5], 5) == 3 // 1
        @test_throws ArgumentError relative_expansion(SimpleGraph(path_graph(25)),
                                                       [1, 2], 1)
        # Single-vertex port subset → vacuous constraint → +∞.
        @test relative_expansion(SimpleGraph(path_graph(4)), [1], 4) ==
              typemax(Int) // 1
        @test relative_expansion(SimpleGraph(complete_graph(5)), [3], 5) ==
              typemax(Int) // 1
    end

    @testset "verify_desideratum_4" begin
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        aux = build_initial_aux_graph(z, as_css(Surface(3, 3)))
        result = verify_desideratum_4(aux, 3)
        @test result isa Bool
        port_vertices = sort(collect(Set(reduce(vcat, aux.port_function))))
        β = relative_expansion(aux.graph, port_vertices, 3)
        @test result == (β ≥ 1)
    end

    @testset "build_initial_aux_graph: disconnected → @warn + auto-repair" begin
        # Hx with two stabs overlapping the support on disjoint qubit pairs.
        Hx = falses(2, 5)
        Hx[1, 1] = Hx[1, 2] = true
        Hx[2, 3] = Hx[2, 4] = true
        Hz = falses(0, 5)
        code = CSS(Matrix{Bool}(Hx), Matrix{Bool}(Hz))
        aux = @test_logs (:warn, r"disconnected"i) match_mode=:any begin
            build_initial_aux_graph([1, 2, 3, 4], code)
        end
        @test is_connected(aux.graph)
        @test nv(aux.graph) == 4
        @test ne(aux.graph) == 3
        @test has_edge(aux.graph, 1, 3)
    end

    @testset "Merge helpers" begin
        for n in (2, 3, 5, 8)
            H = canonical_H_R(n)
            @test size(H) == (n - 1, n)
            for i in 1:(n - 1)
                @test H[i, i] && H[i, i + 1]
                @test sum(H[i, :]) == 2
            end
        end

        g = SimpleGraph(cycle_graph(4))
        st = skiptree(g)
        bridge = adapter_bridge_x_check_matrix(st, st, 4)
        @test size(bridge) == (3, size(st.T, 2) + 4 + size(st.T, 2))
        @test_throws ArgumentError adapter_bridge_x_check_matrix(st, st, 1)
    end

    @testset "port_function_as_matrix" begin
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        c = as_css(Surface(3, 3))
        n_qubits = size(c.Hx, 2)
        aux = build_initial_aux_graph(z, c)
        F = port_function_as_matrix(aux, n_qubits)
        @test size(F) == (nv(aux.graph), n_qubits)
        for q in z
            @test sum(F[:, q]) == 1
        end
        for q in 1:n_qubits
            q ∉ z && @test sum(F[:, q]) == 0
        end
        for (i, q) in enumerate(z)
            @test F[i, q] == true
        end
    end

    @testset "incidence_matrix_transpose" begin
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        aux = build_initial_aux_graph(z, as_css(Surface(3, 3)))
        GT = incidence_matrix_transpose(aux)
        @test size(GT) == (nv(aux.graph), ne(aux.graph))
        @test all(sum(GT; dims = 1) .== 2)
        for (e_idx, edge) in enumerate(edges(aux.graph))
            @test GT[src(edge), e_idx] && GT[dst(edge), e_idx]
        end
    end

    @testset "stabilizer_modification_matrix" begin
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        c = as_css(Surface(3, 3))
        aux = build_initial_aux_graph(z, c)
        n_stabs = size(c.Hx, 1)
        M = stabilizer_modification_matrix(aux, n_stabs)
        @test size(M) == (n_stabs, ne(aux.graph))
        for s in 1:n_stabs
            @test sum(M[s, :]) == length(aux.stabilizer_matchings[s])
        end
        @test sum(M) == sum(length(ms) for ms in aux.stabilizer_matchings)
        @test_throws DimensionMismatch stabilizer_modification_matrix(aux, n_stabs + 1)
    end

    @testset "assemble_merged_code argument validation" begin
        c_s = as_css(Surface(3, 3))
        z_s = sort(findall(!iszero,
                           stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        aux_s = cellulate_long_cycles!(
            build_initial_aux_graph(z_s, c_s);
            max_cycle_len = Int(maximum(sum(c_s.Hx; dims = 2))))
        st_s = skiptree(aux_s.graph)

        c_bb = bb_2bga(6, 6, [(:x, 3), (:y, 1), (:y, 2)],
                              [(:y, 3), (:x, 1), (:x, 2)])
        n_bb = code_n(c_bb)
        c_bb_css = as_css(c_bb)
        z_bb = sort(findall(!iszero,
                            stab_to_gf2(logz_ops(c_bb))[1, n_bb + 1:2n_bb]))
        @assert length(z_s) != length(z_bb)
        aux_bb = cellulate_long_cycles!(
            build_initial_aux_graph(z_bb, c_bb_css);
            max_cycle_len = Int(maximum(sum(c_bb_css.Hx; dims = 2))))
        st_bb = skiptree(aux_bb.graph)

        @test_throws DimensionMismatch assemble_merged_code_intercode(
            c_s, c_bb_css, aux_s, aux_bb, st_s, st_bb)
    end

    @testset "edge_pair_to_index" begin
        g = SimpleGraph(cycle_graph(5))
        for (k, e) in enumerate(edges(g))
            u, v = minmax(src(e), dst(e))
            @test edge_pair_to_index(g, (u, v)) == k
        end
        @test edge_pair_to_index(SimpleGraph(path_graph(4)), (1, 4)) == 0
        @test_throws ArgumentError edge_pair_to_index(g, (3, 1))
    end

    @testset "Surface(3,3) × Surface(3,3) inter-code → [[33, 1, 3]]" begin
        c1 = as_css(Surface(3, 3))
        c2 = as_css(Surface(3, 3))
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        adapter = build_adapter(CodePair(c1, c2, z, z))
        @test code_n(adapter) == 33                            # 13 + 2 + 3 + 2 + 13
        @test code_k(adapter) == 1
        @test css_commutes(adapter.merged_code)
        @test adapter.adapter_width == 3
        @test distance(adapter, DistanceMIPAlgorithm(solver = HiGHS)) == 3
    end

    @testset "BB[[72, 12, 6]] intra-code joint Z measurement" begin
        c_orig = bb_2bga(6, 6, [(:x, 3), (:y, 1), (:y, 2)],
                              [(:y, 3), (:x, 1), (:x, 2)])
        @test (code_n(c_orig), code_k(c_orig)) == (72, 12)
        c = as_css(c_orig)
        lz = stab_to_gf2(logz_ops(c_orig))
        n0 = code_n(c_orig)
        z1 = sort(findall(!iszero, lz[1, n0 + 1:2n0]))
        z2 = sort(findall(!iszero, lz[2, n0 + 1:2n0]))
        adapter = build_adapter(CodePair(c, c, z1, z2))
        @test css_commutes(adapter.merged_code)
        @test code_k(adapter) == code_k(c_orig) - 1
        # paper Theorem 5 / Table III weight bounds
        @test maximum(sum(parity_matrix_x(adapter), dims = 2)) ≤ 8
        @test maximum(sum(parity_matrix_z(adapter), dims = 2)) ≤ 8
    end

    @testset "Intra-code adapter with unequal-width logicals (w_l ≠ w_r)" begin
        # Trigger the `assemble_merged_code_intracode` T/P trimming path.
        bb = BivariateBicycleViaCirculantMat(12, 6,
                [(:x, 3), (:y, 1), (:y, 2)],
                [(:y, 3), (:x, 1), (:x, 2)])
        n0 = code_n(bb)
        lz = stab_to_gf2(logz_ops(bb))
        weights = [count(!iszero, lz[i, n0 + 1:2n0]) for i in axes(lz, 1)]
        distinct = sort(unique(weights))
        @assert length(distinct) ≥ 2
        i_short = findfirst(==(distinct[1]),    weights)
        i_long  = findfirst(==(distinct[end]),  weights)
        z_short = sort(findall(!iszero, lz[i_short, n0 + 1:2n0]))
        z_long  = sort(findall(!iszero, lz[i_long,  n0 + 1:2n0]))
        @test length(z_short) < length(z_long)

        c = as_css(bb)
        adapter = build_adapter(CodePair(c, c, z_short, z_long))

        @test adapter.adapter_width == length(z_short)
        @test nv(adapter.aux_l.graph) == length(z_short)
        @test nv(adapter.aux_r.graph) == length(z_long)
        @test code_n(adapter) == n0 + ne(adapter.aux_l.graph) +
                                      adapter.adapter_width +
                                      ne(adapter.aux_r.graph)
        @test css_commutes(adapter.merged_code)
        @test code_k(adapter) == code_k(bb) - 1
    end

    @testset "Gross [[144, 12, 12]] BB structural regression" begin
        # Bravyi-et-al-2024 gross BB. Distance MIP at d=12 is out-of-band.
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

        expected_n = n0 + ne(adapter.aux_l.graph) +
                          adapter.adapter_width +
                          ne(adapter.aux_r.graph)
        @test code_n(adapter) == expected_n
        @test adapter.adapter_width == min(length(z1), length(z2))
        @test css_commutes(adapter.merged_code)
        @test code_k(adapter) == code_k(g_orig) - 1
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
        aux = cellulate_long_cycles!(build_initial_aux_graph(z, c);
                                      max_cycle_len = Int(maximum(sum(c.Hx; dims = 2))))
        @test code_n(dc) == code_n(c) + ne(aux.graph)
        @test code_k(dc) == 0
        @test css_commutes(dc)
    end

    @testset "build_adapter dispatcher" begin
        c  = as_css(Surface(3, 3))
        c1 = as_css(Surface(3, 3))
        c2 = as_css(Surface(3, 3))
        z = sort(findall(!iszero,
                         stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        @test code_n(build_adapter(CodePair(c1, c2, z, z))) == 33
        @test build_adapter(CodePair(c, c, z, z)).adapter_width == 3
        @test_throws ArgumentError build_adapter_intercode(CodePair(c, c, z, z))
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
        @test distance(Steane7(), DistanceMIPAlgorithm(solver = HiGHS)) == 3
        @test distance(Surface(3, 3), DistanceMIPAlgorithm(solver = HiGHS)) == 3
    end
end
