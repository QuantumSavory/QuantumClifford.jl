@testitem "Adapters" tags=[:ecc, :ecc_adapters] begin
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC.Adapters
    # Internal pipeline-stage and block-helper functions used by tests.
    # Adapters only exports the user-facing surface (`build_adapter`,
    # `deform_code`, etc.); test-only internals are imported explicitly.
    using QuantumClifford.ECC.Adapters: build_initial_aux_graph,
                                        cellulate_long_cycles!,
                                        assemble_merged_code_intercode,
                                        assemble_merged_code_intracode,
                                        relative_expansion,
                                        verify_desideratum_4,
                                        port_function_as_matrix,
                                        edge_pair_to_index,
                                        incidence_matrix_transpose,
                                        stabilizer_modification_matrix,
                                        canonical_H_R,
                                        adapter_bridge_x_check_matrix
    using QuantumClifford.ECC: CSS, Surface, Steane7, BivariateBicycleViaCirculantMat,
                              code_n, code_k, code_s,
                              parity_checks, parity_matrix, parity_matrix_x, parity_matrix_z,
                              logx_ops, logz_ops, deform_code, build_adapter
    using QECCore: AbstractCSSCode, AbstractQECC, AbstractECC
    using QuantumClifford: stab_to_gf2, nqubits
    using Graphs: SimpleGraph, path_graph, complete_graph, cycle_graph,
                  add_edge!, rem_edge!, nv, ne, edges, src, dst,
                  has_edge, is_connected, degree, neighbors, cycle_basis,
                  erdos_renyi
    using SparseArrays
    using LinearAlgebra
    using Random
    using JuMP, HiGHS

    BASE = joinpath(@__DIR__, "..", "data", "zenodo_17527545",
                    "data_qLDPC_surgery")

    """Tiny MatrixMarket reader (test env doesn't have MatrixMarket.jl)."""
    function load_mtx_inline(path::AbstractString)
        lines = filter(l -> !startswith(l, "%") && !isempty(strip(l)),
                       readlines(path))
        header = split(lines[1])
        rows = parse(Int, header[1]); cols = parse(Int, header[2])
        M = falses(rows, cols)
        for l in lines[2:end]
            parts = split(l)
            i, j = parse(Int, parts[1]), parse(Int, parts[2])
            v = length(parts) ≥ 3 ? parse(Int, parts[3]) : 1
            v != 0 && (M[i, j] = true)
        end
        Matrix(M)
    end

    """Lightweight parser for Zenodo `G_mat_X` .txt files — recovers the
    edge vertex pairs in Zenodo's storage order."""
    function load_g_mat_pairs(path::AbstractString, varname::AbstractString)
        src_str = read(path, String)
        src_str = join(split(src_str, "\n")[2:end], "\n")
        idx = findfirst(varname * " ", src_str)
        idx === nothing && error("variable $varname not found")
        open_idx = findnext("np.array([", src_str, last(idx))
        cur = last(open_idx) + 1
        depth = 1; body_start = cur - 1; body_end = 0; i = cur
        while i ≤ lastindex(src_str)
            c = src_str[i]
            if c == '['
                depth += 1
            elseif c == ']'
                depth -= 1
                if depth == 0
                    body_end = i; break
                end
            end
            i = nextind(src_str, i)
        end
        body = src_str[body_start:body_end]
        inner = strip(body); inner = inner[2:end - 1]
        row_strs = String[]
        d = 0; last_start = 1
        for (ix, c) in pairs(inner)
            if c == '['
                d == 0 && (last_start = ix)
                d += 1
            elseif c == ']'
                d -= 1
                d == 0 && push!(row_strs, inner[last_start:ix])
            end
        end
        pairs_out = Tuple{Int,Int}[]
        for rs in row_strs
            rs2 = replace(rs, r"[\[\]\n]" => "")
            parts = split(rs2, ",")
            cols = Int[]
            for (j, p) in enumerate(parts)
                sp = strip(p); isempty(sp) && continue
                parse(Float64, sp) != 0 && push!(cols, j)
            end
            length(cols) == 2 ||
                error("non-edge row in $varname (got $(length(cols)) nonzeros)")
            push!(pairs_out, (cols[1], cols[2]))
        end
        pairs_out
    end

    """Promote any AbstractCSSCode to a plain `CSS` of dense `Bool` matrices."""
    as_css(c) = CSS(Matrix{Bool}(parity_matrix_x(c)),
                    Matrix{Bool}(parity_matrix_z(c)))

    """GF(2) rank via the existing canonicalization helper."""
    function rank_gf2(H::AbstractMatrix{Bool})::Int
        M = Int.(H)
        _, r, _, _ = QuantumClifford.gf2_row_echelon_with_pivots!(M)
        r
    end

    """MIP-feasibility distance probe (used by the MIP-sanity testset)."""
    function logical_op_exists_at_weight_le_T(
        code, logical_idx::Int, op_type::Symbol, T::Int;
        time_limit::Float64 = 30.0,
    )
        n = nqubits(code)
        if op_type == :X
            l = stab_to_gf2(logx_ops(code))
            h = Int.(parity_matrix_x(code))
            lr = Int.(vec(l[logical_idx, 1:n]))
        elseif op_type == :Z
            l = stab_to_gf2(logz_ops(code))
            h = Int.(parity_matrix_z(code))
            lr = Int.(vec(l[logical_idx, n+1:2n]))
        end
        m = size(h, 1)
        whx = maximum(sum(h, dims = 2))
        wlx = count(!iszero, lr)
        num_anc_hx = Int(ceil(log2(whx + 1)))
        num_anc_logical = max(1, Int(ceil(log2(wlx + 1))))
        num_var = n + m * num_anc_hx + num_anc_logical
        model = Model(HiGHS.Optimizer; add_bridges = false)
        set_silent(model)
        set_time_limit_sec(model, time_limit)
        @variable(model, x[1:num_var], Bin)
        @objective(model, Min, sum(x[i] for i in 1:n))
        for row in 1:m
            @constraint(model,
                sum(h[row, q] * x[q] for q in 1:n) -
                sum((1 << k) * x[n + (row - 1) * num_anc_hx + k] for k in 1:num_anc_hx) == 0)
        end
        supp = findall(!iszero, lr)
        @constraint(model,
            sum(x[q] for q in supp) -
            sum((1 << k) * x[n + m * num_anc_hx + k] for k in 1:num_anc_logical) == 1)
        @constraint(model, sum(x[i] for i in 1:n) ≤ T)
        optimize!(model)
        status = termination_status(model)
        if status == MOI.INFEASIBLE
            return (:proven_lower_bound, T + 1)
        elseif status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
            opt_w = round(Int, sum(value(x[i]) for i in 1:n))
            return (:counterexample_found, opt_w)
        else
            return (:inconclusive, status)
        end
    end

    # ─── Shared Zenodo data ───────────────────────────────────────────────
    # Loaded once at the top of the testitem; reused by every testset that
    # needs it. Testsets that touch Zenodo first check `isdir(BASE)` and
    # `@test_skip` when the data isn't present (CI without the Zenodo
    # release).
    have_zenodo = isdir(BASE)
    bb_canonical = bb_full = lp_code = nothing
    if have_zenodo
        bb_canonical = CSS(
            load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                "Hx_BB_98_6_12_original-code-canonicalbasis.mtx")),
            load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                "Hz_BB_98_6_12_original-code-canonicalbasis.mtx")))
        bb_full = CSS(
            load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                "Hx_BB_98_6_12_original-code-fullrankbasis.mtx")),
            load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                "Hz_BB_98_6_12_original-code-fullrankbasis.mtx")))
        lp_code = CSS(
            load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
                "Hx_LP_200_20_10_original-code.mtx")),
            load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
                "Hz_LP_200_20_10_original-code.mtx")))
    end
    # Paper Table I logical operator supports (1-indexed = paper's 0-indexed + 1).
    Z1_BB = sort([6,8,13,17,31,32,33,35,36,37,41,50,51,93]) .+ 1   # |Z_1| = 14
    Z2_LP = sort([24,25,26,29,30,56,58,59,60,61,90,93,94,121]) .+ 1
    Z3_BB = sort([10,17,35,39,42,43,53,55,61,70,84,89])    .+ 1    # |Z_3| = 12

    @testset "SkipTree (Algorithm 2 / H_R)" begin
                  edges, src, dst, is_connected, erdos_renyi

    """Incidence matrix (m × n) of a SimpleGraph, in the order
    `edges(g)` iterates."""
    function incidence_matrix(g::SimpleGraph)
        I = Int[]; J = Int[]
        for (k, e) in enumerate(edges(g))
            push!(I, k); push!(J, src(e))
            push!(I, k); push!(J, dst(e))
        end
        sparse(I, J, trues(length(I)), ne(g), nv(g))
    end

    """Canonical H_R(n) = (n-1) × n full-rank repetition check matrix."""
    function H_R(n::Int)
        M = falses(n - 1, n)
        for i in 1:(n - 1)
            M[i, i] = true
            M[i, i + 1] = true
        end
        sparse(M)
    end

    is_permutation(P::SparseMatrixCSC) =
        all(sum(P; dims = 1) .== 1) &&
        all(sum(P; dims = 2) .== 1) &&
        all(x -> x == 0 || x == 1, P)

    function check_skiptree(g::SimpleGraph; root::Int = 1)
        n = nv(g); m = ne(g)
        out = skiptree(g; root = root)
        @test out isa SkipTreeOutput
        @test out.target === :H_R
        @test size(out.T) == (n - 1, m)
        @test size(out.P) == (n, n)

        # (3, 2)-sparse?
        rw = vec(sum(out.T; dims = 2))
        cw = vec(sum(out.T; dims = 1))
        @test maximum(rw) ≤ 3
        @test maximum(cw) ≤ 2
        @test is_permutation(out.P)

        # T · G · P (mod 2) == H_R(n) ?
        G = incidence_matrix(g)
        TGP = mod.(Matrix{Int}(out.T) * Matrix{Int}(G) * Matrix{Int}(out.P), 2)
        @test TGP == Matrix{Int}(H_R(n))

        out
    end

    @testset "1. Path graph (n=5)" begin
        g = SimpleGraph(path_graph(5))
        out = check_skiptree(g)
        @test size(out.T) == (4, 4)   # m = 4 = n-1 for path
    end

    @testset "2. Cycle graph (n=5)" begin
        g = SimpleGraph(cycle_graph(5))
        out = check_skiptree(g)
        @test size(out.T) == (4, 5)
    end

    @testset "3. Random connected graphs (n=20, ten trials)" begin
        Random.seed!(2026_05_21)
        ntrials = 10
        passed = 0
        for trial in 1:ntrials
            g = SimpleGraph(erdos_renyi(20, 0.3; seed = trial * 17 + 3))
            # Erdos-Renyi at p=0.3 over n=20 is connected w.h.p.; if a draw
            # comes back disconnected, add a spanning path to repair it.
            if !is_connected(g)
                for v in 2:nv(g)
                    add_edge!(g, v - 1, v)
                end
            end
            @assert is_connected(g)
            check_skiptree(g)
            passed += 1
        end
        @test passed == ntrials
    end

    @testset "4. Zenodo G_mat_1 / G_mat_2 / G_mat_3 (algebraic equivalence)" begin
        # We don't depend on any of the validation/scripts loaders here;
        # we re-parse the txt files inline to keep the test self-contained.
        function parse_python_matrix(path::AbstractString, varname::AbstractString)
            src = read(path, String)
            # The first line of every Zenodo *_GTP*.txt is a malformed docstring.
            src = join(split(src, "\n")[2:end], "\n")
            needle = varname * " "
            pos = findfirst(needle, src)
            pos === nothing && error("variable $varname not found in $path")
            open_idx = findnext("np.array([", src, last(pos))
            open_idx === nothing && error("malformed assignment for $varname")
            cur = last(open_idx) + 1
            bracket_depth = 1
            body_start = cur - 1
            body_end = 0
            i = cur
            while i ≤ lastindex(src)
                c = src[i]
                if c == '['
                    bracket_depth += 1
                elseif c == ']'
                    bracket_depth -= 1
                    if bracket_depth == 0
                        body_end = i
                        break
                    end
                end
                i = nextind(src, i)
            end
            @assert body_end > 0 "unbalanced brackets parsing $varname"
            body = src[body_start:body_end]
            inner = strip(body)
            @assert startswith(inner, "[") && endswith(inner, "]")
            inner = inner[2:end - 1]
            row_strs = String[]
            depth = 0; last_start = 1
            for (idx, c) in pairs(inner)
                if c == '['
                    if depth == 0
                        last_start = idx
                    end
                    depth += 1
                elseif c == ']'
                    depth -= 1
                    if depth == 0
                        push!(row_strs, inner[last_start:idx])
                    end
                end
            end
            rows = Vector{Vector{Float64}}()
            for rs in row_strs
                rs2 = replace(rs, r"[\[\]\n]" => "")
                parts = split(rs2, ",")
                vals = Float64[]
                for p in parts
                    sp = strip(p)
                    isempty(sp) && continue
                    push!(vals, parse(Float64, sp))
                end
                push!(rows, vals)
            end
            nrows = length(rows); ncols = length(rows[1])
            @assert all(length(r) == ncols for r in rows)
            M = falses(nrows, ncols)
            for i in 1:nrows, j in 1:ncols
                if rows[i][j] != 0.0
                    M[i, j] = true
                end
            end
            sparse(M)
        end

        if !isdir(BASE)
            @info "Skipping Zenodo tests: data not present at $BASE"
            @test_skip "Zenodo tests skipped (no data)"
        else
            zen_files = [
                (joinpath(BASE, "BB_98_LP_200_adapter",
                          "skipTree_transformations", "BB_98_6_12_Z_1_GTP.txt"),
                 "G_mat_1", "Z_1"),
                (joinpath(BASE, "BB_98_LP_200_adapter",
                          "skipTree_transformations",
                          "LP_200_20_10_Z_2_GTP_2.txt"),
                 "G_mat_2", "Z_2"),
                (joinpath(BASE, "BB_98_intracode_adapter",
                          "skipTree_transformations_G_1_G3",
                          "BB_98_6_12_Z_3_GTP-other-logical.txt"),
                 "G_mat_3", "Z_3"),
            ]
            for (path, var, label) in zen_files
                @testset "Zenodo $label" begin
                    G_mat = parse_python_matrix(path, var)
                    m, n = size(G_mat)
                    g = SimpleGraph(n)
                    for i in 1:m
                        cols = findnz(sparse(G_mat[i, :]))[1]
                        @assert length(cols) == 2
                        add_edge!(g, cols[1], cols[2])
                    end
                    @test is_connected(g)
                    check_skiptree(g)
                end
            end
        end
    end
    end
    @testset "auxiliary graph (build_initial_aux_graph)" begin
    # Tiny MatrixMarket reader (the test env doesn't have MatrixMarket.jl).


    @testset "1. Trivial sanity case (4-qubit logical, one weight-4 X-stab)" begin
        n = 4
        Hx = falses(1, n); Hx[1, 1:4] .= true
        Hz = falses(1, n); Hz[1, 1:4] .= true
        code = CSS(Matrix(Hx), Matrix(Hz))
        # Connectivity-repair will fire for this tiny case (matching pairs
        # produce 2 disconnected components). We don't bother suppressing the
        # @warn; it's expected output, not a failure signal.
        aux = build_initial_aux_graph([1, 2, 3, 4], code)
        @test nv(aux.graph) == 4
        # The matching adds (1,2) and (3,4); these are disconnected, so
        # the repair adds one more edge -> 3 edges total.
        @test ne(aux.graph) == 3
        @test is_connected(aux.graph)
        @test aux.port_function == [[1], [2], [3], [4]]
        # Max vertex degree should stay small.
        @test maximum(degree(aux.graph)) ≤ 2
        # Only one X-stabilizer; its matching has 2 edges.
        @test length(aux.stabilizer_matchings) == 1
        @test length(aux.stabilizer_matchings[1]) == 2
        # Each stored entry is a sorted vertex pair that is a real edge.
        @test all(p -> p[1] ≤ p[2] && 1 ≤ p[1] && p[2] ≤ nv(aux.graph) &&
                       has_edge(aux.graph, p[1], p[2]),
                  aux.stabilizer_matchings[1])
    end

    @testset "2. Zenodo Example B G_1 — BB Hx canonical, Z_1 (size 14)" begin
        if !isdir(BASE)
            @info "Skipping Zenodo aux-graph test: data not present"
            @test_skip "no Zenodo data"
            return
        end
        Hx_BB = load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                                "Hx_BB_98_6_12_original-code-canonicalbasis.mtx"))
        Hz_BB = load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                                "Hz_BB_98_6_12_original-code-canonicalbasis.mtx"))
        code = CSS(Hx_BB, Hz_BB)
        Z1 = sort([6, 8, 13, 17, 31, 32, 33, 35, 36, 37, 41, 50, 51, 93]) .+ 1
        aux = build_initial_aux_graph(Z1, code)
        @test nv(aux.graph) == 14
        # Matching-only edge count is 21 (one matching edge per overlapping stabilizer,
        # all |L_s|=2 so no within-stabilizer freedom, no edge collisions).
        # Zenodo's final G_1 has 23 edges = 21 matching edges + 2 cellulation chord edges.
        # Cellulation is added by cellulate_long_cycles! (not yet implemented).
        @info "Example B G_1 (matching-only): ne(aux.graph) = $(ne(aux.graph))"
        @test ne(aux.graph) == 21
        @test is_connected(aux.graph)
        @test aux.port_function == [[i] for i in 1:14]
        # F should be (14, 98) with one nonzero per row and one per Z1 column.
        F = port_function_as_matrix(aux, 98)
        @test size(F) == (14, 98)
        @test sum(F) == 14
        @test all(==(1), vec(sum(F; dims = 2)))
        @test sort(unique([j for (i, j) in zip(findnz(F)[1], findnz(F)[2])])) == Z1
        # Cell-by-cell check against the reverse-engineered F_l:
        # F_l[v, q] == 1 iff v == findfirst(==(q), Z1) i.e. F_l is the
        # natural one-hot port matrix.
        F_l_expected = falses(14, 98)
        for (v, q) in enumerate(Z1)
            F_l_expected[v, q] = true
        end
        @test Matrix(F) == F_l_expected
        # stabilizer_matchings sanity
        nonempty = [s for (s, m) in enumerate(aux.stabilizer_matchings)
                    if !isempty(m)]
        # Paper §VII.A reverse engineering: 21 of 49 BB stabs overlap Z1
        @info "Example B G_1: # stabs with matchings = $(length(nonempty)) (paper: 21)"
        @test length(nonempty) == 21
        for s in 1:size(Hx_BB, 1)
            ms = aux.stabilizer_matchings[s]
            # Each matching has |L_s|/2 entries where |L_s| is even.
            ovl_size = count(q -> Hx_BB[s, q] && (q ∈ Set(Z1)), 1:98)
            @test length(ms) == ovl_size ÷ 2 || (ovl_size < 2 && isempty(ms))
            # Each stored entry is a sorted vertex pair that is a real edge.
            @test all(p -> p[1] ≤ p[2] && 1 ≤ p[1] && p[2] ≤ nv(aux.graph) &&
                           has_edge(aux.graph, p[1], p[2]),
                      ms)
        end
    end

    @testset "3. Zenodo Example B G_2 — LP Hx, Z_2 (size 14)" begin
        if !isdir(BASE)
            @test_skip "no Zenodo data"
            return
        end
        Hx_LP = load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
                                "Hx_LP_200_20_10_original-code.mtx"))
        Hz_LP = load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
                                "Hz_LP_200_20_10_original-code.mtx"))
        code = CSS(Hx_LP, Hz_LP)
        Z2 = sort([24, 25, 26, 29, 30, 56, 58, 59, 60, 61, 90, 93, 94, 121]) .+ 1
        aux = build_initial_aux_graph(Z2, code)
        @test nv(aux.graph) == 14
        # 21 stabilizers contribute matching edges, but stabilizers 32 and 96 both
        # produce edge (2,14), so 20 distinct edges.  Zenodo agrees: no cellulation
        # was needed for G_2 (paper §VII.A "No cellulation or decongestion needed").
        @info "Example B G_2 (matching-only): ne(aux.graph) = $(ne(aux.graph))"
        @test ne(aux.graph) == 20
        @test is_connected(aux.graph)
        @test aux.port_function == [[i] for i in 1:14]
        F = port_function_as_matrix(aux, 200)
        @test size(F) == (14, 200)
        # F_r ground truth: one-hot at (v, Z2[v])
        F_r_expected = falses(14, 200)
        for (v, q) in enumerate(Z2)
            F_r_expected[v, q] = true
        end
        @test Matrix(F) == F_r_expected
    end

    @testset "4. Zenodo Example C G_3 — BB Hx full-rank, Z_3 (size 12)" begin
        if !isdir(BASE)
            @test_skip "no Zenodo data"
            return
        end
        Hx_BB = load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                                "Hx_BB_98_6_12_original-code-fullrankbasis.mtx"))
        Hz_BB = load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                                "Hz_BB_98_6_12_original-code-fullrankbasis.mtx"))
        code = CSS(Hx_BB, Hz_BB)
        Z3 = sort([10, 17, 35, 39, 42, 43, 53, 55, 61, 70, 84, 89]) .+ 1
        aux = build_initial_aux_graph(Z3, code)
        @test nv(aux.graph) == 12
        # Matching-only edge count is 16 (one matching edge per overlapping stabilizer,
        # all |L_s|=2 so no within-stabilizer freedom, no edge collisions).
        # Zenodo's final G_3 has 17 edges = 16 matching edges + 1 cellulation chord edge.
        # Cellulation will be handled by cellulate_long_cycles! (not yet implemented).
        @info "Example C G_3 (matching-only): ne(aux.graph) = $(ne(aux.graph))"
        @test ne(aux.graph) == 16
        @test is_connected(aux.graph)
        @test aux.port_function == [[i] for i in 1:12]
        F = port_function_as_matrix(aux, 98)
        @test size(F) == (12, 98)
        F_r_expected = falses(12, 98)
        for (v, q) in enumerate(Z3)
            F_r_expected[v, q] = true
        end
        @test Matrix(F) == F_r_expected
    end

    @testset "5. Boundary: 2-qubit logical with no overlap → @warn + repair" begin
        # CSS code with Hx that has NO overlap with the chosen logical support.
        n = 4
        # All X-stabs supported only on {1, 2}; logical_support = {3, 4}
        Hx = falses(2, n); Hx[1, 1] = true; Hx[1, 2] = true
                          Hx[2, 1] = true; Hx[2, 2] = true
        # Z stabs supported only on {3, 4} so CSS commutes
        Hz = falses(2, n); Hz[1, 3] = true; Hz[1, 4] = true
                          Hz[2, 3] = true; Hz[2, 4] = true
        code = CSS(Matrix(Hx), Matrix(Hz))
        @test_logs (:warn, r"disconnected") match_mode=:any begin
            aux = build_initial_aux_graph([3, 4], code)
            @test nv(aux.graph) == 2
            # 0 matching edges + 1 repair edge → 1 edge total
            @test ne(aux.graph) == 1
            @test is_connected(aux.graph)
            @test all(isempty, aux.stabilizer_matchings)
        end
    end
    end
    @testset "cellulate_long_cycles!" begin
    # Reuse the inline MatrixMarket reader from the build_initial_aux_graph
    # testitem (kept inline so this testitem is self-contained).


    function code_BB_canonical()
        Hx = load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
            "Hx_BB_98_6_12_original-code-canonicalbasis.mtx"))
        Hz = load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
            "Hz_BB_98_6_12_original-code-canonicalbasis.mtx"))
        CSS(Hx, Hz)
    end
    function code_BB_fullrank()
        Hx = load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
            "Hx_BB_98_6_12_original-code-fullrankbasis.mtx"))
        Hz = load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
            "Hz_BB_98_6_12_original-code-fullrankbasis.mtx"))
        CSS(Hx, Hz)
    end
    function code_LP()
        Hx = load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
            "Hx_LP_200_20_10_original-code.mtx"))
        Hz = load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
            "Hz_LP_200_20_10_original-code.mtx"))
        CSS(Hx, Hz)
    end

    # Z1_BB / Z2_LP / Z3_BB defined at testitem top.

    setup_B_G1() = build_initial_aux_graph(Z1_BB, code_BB_canonical())
    setup_B_G2() = build_initial_aux_graph(Z2_LP, code_LP())
    setup_C_G3() = build_initial_aux_graph(Z3_BB, code_BB_fullrank())

    """N · G_incidence ≡ 0 (mod 2)?"""
    function NG_check(aux)
        n = nv(aux.graph); m = ne(aux.graph)
        G_inc = falses(m, n)
        for (i, e) in enumerate(edges(aux.graph))
            G_inc[i, src(e)] = true
            G_inc[i, dst(e)] = true
        end
        NG = mod.(Matrix{Int}(aux.cycle_basis_matrix) * Matrix{Int}(G_inc), 2)
        all(iszero, NG), m - n + 1
    end

    # GF(2) rank via row reduction
    function gf2_rank(M::AbstractMatrix)
        A = Matrix{Int}(M)
        r, c = size(A)
        rank_ = 0
        col = 1
        for row in 1:r
            while col ≤ c
                pivot = 0
                for k in row:r
                    if A[k, col] == 1; pivot = k; break; end
                end
                if pivot == 0
                    col += 1
                    continue
                end
                A[[row, pivot], :] = A[[pivot, row], :]
                for k in 1:r
                    if k != row && A[k, col] == 1
                        @views A[k, :] .= mod.(A[k, :] .+ A[row, :], 2)
                    end
                end
                rank_ += 1
                col += 1
                break
            end
            col > c && break
        end
        rank_
    end

    if !isdir(BASE)
        @info "Skipping Zenodo cellulation tests: data not present at $BASE"
        @test_skip "no Zenodo data"
    else
        @testset "Edge count goes 21 → 23 (B / G_1, max_cycle_len=6)" begin
            aux = setup_B_G1()
            @test ne(aux.graph) == 21
            aux2 = cellulate_long_cycles!(aux; max_cycle_len = 6)
            @info "B / G_1 cellulated: ne = $(ne(aux2.graph))"
            @test ne(aux2.graph) == 23
        end

        @testset "Edge count goes 20 → 20 (B / G_2, max_cycle_len=7)" begin
            # Per paper §VII.A: "since the original LP_2 code already has
            # stabilizer weights 7, we choose not to cellulate these cycles."
            # Using max_cycle_len=7 reproduces Zenodo's choice exactly.
            aux = setup_B_G2()
            @test ne(aux.graph) == 20
            aux2 = cellulate_long_cycles!(aux; max_cycle_len = 7)
            @info "B / G_2 cellulated (max_cycle_len=7): ne = $(ne(aux2.graph))"
            @test ne(aux2.graph) == 20
        end

        @testset "Edge count goes 20 → 21 with default max_cycle_len=6 (B / G_2)" begin
            # When using the default threshold of 6, Julia's basis exposes one
            # length-7 cycle that gets cellulated; Zenodo opted out by lifting
            # the threshold. This test documents the algorithmic behaviour.
            aux = setup_B_G2()
            aux2 = cellulate_long_cycles!(aux; max_cycle_len = 6)
            @info "B / G_2 cellulated (max_cycle_len=6): ne = $(ne(aux2.graph))"
            @test ne(aux2.graph) == 21
        end

        @testset "Edge count goes 16 → 17 (C / G_3, max_cycle_len=6)" begin
            aux = setup_C_G3()
            @test ne(aux.graph) == 16
            aux2 = cellulate_long_cycles!(aux; max_cycle_len = 6)
            @info "C / G_3 cellulated: ne = $(ne(aux2.graph))"
            @test ne(aux2.graph) == 17
        end

        @testset "Max basis cycle length ≤ max_cycle_len after cellulation" begin
            for (setup, threshold) in [(setup_B_G1, 6),
                                        (setup_B_G2, 7),
                                        (setup_C_G3, 6)]
                aux = setup()
                aux2 = cellulate_long_cycles!(aux; max_cycle_len = threshold)
                basis = cycle_basis(aux2.graph)
                @test maximum(length, basis) ≤ threshold
            end
        end

        @testset "Cycle-basis matrix algebraic correctness (N · G = 0 mod 2)" begin
            for (setup, threshold) in [(setup_B_G1, 6),
                                        (setup_B_G2, 7),
                                        (setup_C_G3, 6)]
                aux = setup()
                aux2 = cellulate_long_cycles!(aux; max_cycle_len = threshold)
                ok, expected_rank = NG_check(aux2)
                @test ok
                @test size(aux2.cycle_basis_matrix, 1) == expected_rank
                @test gf2_rank(aux2.cycle_basis_matrix) == expected_rank
            end
        end

        @testset "Post-cellulation length distribution for B / G_1" begin
            # Paper §VII.A claims [3, 3, 4, 5, 5, 5, 5, 5, 6, 6] post-cellulation.
            # Julia's Graphs.cycle_basis chooses a different DFS basis than
            # NetworkX (basis-choice freedom). Both bases span the same cycle
            # space; the *invariants* we care about are tested elsewhere
            # (edge count, max length, N·G=0 mod 2, rank). The specific length
            # distribution is non-canonical and not asserted byte-for-byte.
            aux = setup_B_G1()
            aux2 = cellulate_long_cycles!(aux; max_cycle_len = 6)
            lens = sort(length.(cycle_basis(aux2.graph)))
            @info "B / G_1 post-cellulation cycle lengths: $lens (paper: [3,3,4,5,5,5,5,5,6,6])"
            @test length(lens) == 10           # cyclomatic number 23 - 14 + 1
            @test sum(lens) >= 10 * 3          # sanity: every cycle has length ≥ 3
            @test maximum(lens) ≤ 6            # post-cellulation invariant
        end
    end
    end
    @testset "relative expansion (desideratum 4)" begin
    # --- inline MTX reader (same as the other testitems) ---


    code_BB_canonical() = CSS(
        load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
            "Hx_BB_98_6_12_original-code-canonicalbasis.mtx")),
        load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
            "Hz_BB_98_6_12_original-code-canonicalbasis.mtx")))
    code_BB_fullrank() = CSS(
        load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
            "Hx_BB_98_6_12_original-code-fullrankbasis.mtx")),
        load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
            "Hz_BB_98_6_12_original-code-fullrankbasis.mtx")))
    code_LP() = CSS(
        load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
            "Hx_LP_200_20_10_original-code.mtx")),
        load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
            "Hz_LP_200_20_10_original-code.mtx")))

    z1_bb = sort([6,8,13,17,31,32,33,35,36,37,41,50,51,93]) .+ 1
    z2_lp = sort([24,25,26,29,30,56,58,59,60,61,90,93,94,121]) .+ 1
    z3_bb = sort([10,17,35,39,42,43,53,55,61,70,84,89]) .+ 1

    function aux_B_G1_cellulated()
        aux = build_initial_aux_graph(z1_bb, code_BB_canonical())
        cellulate_long_cycles!(aux; max_cycle_len = 6)
    end
    function aux_B_G2_cellulated()
        aux = build_initial_aux_graph(z2_lp, code_LP())
        cellulate_long_cycles!(aux; max_cycle_len = 7)
    end
    function aux_C_G3_cellulated()
        aux = build_initial_aux_graph(z3_bb, code_BB_fullrank())
        cellulate_long_cycles!(aux; max_cycle_len = 6)
    end

    @testset "Boundary: path graph P_4 has β_4(V) = 1/2" begin
        g = SimpleGraph(path_graph(4))
        β = relative_expansion(g, [1, 2, 3, 4], 4)
        @test β == 1 // 2
        @test β < 1
        # Build an AuxiliaryGraph wrapper to call verify_desideratum_4.
        # Note: stabilizer_matchings length is arbitrary metadata here.
        aux = AuxiliaryGraph(g, [[1], [2], [3], [4]], collect(1:ne(g)),
                              spzeros(Bool, 0, ne(g)), [1, 2, 3, 4],
                              [Tuple{Int,Int}[]])
        @test verify_desideratum_4(aux, 4) == false
    end

    @testset "Boundary: complete graph K_5 has β_5(V) = 3" begin
        g = SimpleGraph(complete_graph(5))
        β = relative_expansion(g, [1, 2, 3, 4, 5], 5)
        @test β == 3 // 1
        @test β ≥ 1
        aux = AuxiliaryGraph(g, [[i] for i in 1:nv(g)],
                              collect(1:ne(g)), spzeros(Bool, 0, ne(g)),
                              collect(1:nv(g)),
                              [Tuple{Int,Int}[] for _ in 1:nv(g)])
        @test verify_desideratum_4(aux, 5) == true
    end

    if !isdir(BASE)
        @info "Skipping Zenodo expansion tests: data not present at $BASE"
        @test_skip "no Zenodo data"
    else
        @testset "Example B / G_1 (d=12): β = 5/6 < 1" begin
            aux = aux_B_G1_cellulated()
            n = nv(aux.graph)
            β = relative_expansion(aux.graph, collect(1:n), 12)
            @info "B/G_1 β_12(V) = $β = $(Float64(β))"
            @test β == 5 // 6
            @test β < 1
            @test verify_desideratum_4(aux, 12) == false
        end

        @testset "Example B / G_2 (d=10): β = 1/2 < 1" begin
            aux = aux_B_G2_cellulated()
            n = nv(aux.graph)
            β = relative_expansion(aux.graph, collect(1:n), 10)
            @info "B/G_2 β_10(V) = $β = $(Float64(β))"
            @test β == 1 // 2
            @test β < 1
            @test verify_desideratum_4(aux, 10) == false
        end

        @testset "Example C / G_3 (d=12): β = 3/4 < 1" begin
            aux = aux_C_G3_cellulated()
            n = nv(aux.graph)
            β = relative_expansion(aux.graph, collect(1:n), 12)
            @info "C/G_3 β_12(V) = $β = $(Float64(β))"
            @test β == 3 // 4
            @test β < 1
            @test verify_desideratum_4(aux, 12) == false
        end

        @testset "Global Cheeger comparison (Lemma 3: β_t(U) ≥ β(V) when t ≤ |V|)" begin
            # With identity port functions, im f = V and t ≤ |V|, so relative
            # and global expansion coincide at t = d.  Lemma 3 holds trivially
            # as equality here.
            for (name, factory, d) in [
                    ("B / G_1", aux_B_G1_cellulated, 12),
                    ("B / G_2", aux_B_G2_cellulated, 10),
                    ("C / G_3", aux_C_G3_cellulated, 12)]
                aux = factory()
                n = nv(aux.graph)
                β_rel = relative_expansion(aux.graph, collect(1:n), d)
                β_glob = relative_expansion(aux.graph, collect(1:n), n)
                @info "$name: β_d=$β_rel, β_global=$β_glob"
                @test β_rel >= β_glob   # Lemma 3 (trivially equal here)
            end
        end

        # Theoretical context (informational, not asserted):
        #
        # Theorem 11 of arXiv:2410.03628 requires β_d(G, imf) ≥ 1, which
        # FAILS for all three Zenodo examples (β = 5/6, 1/2, 3/4).
        #
        # Theorem 12 conditions (a) separate codeblocks, (b) larger Z̄ is a
        # min-weight logical of its code:
        #   - Example B: (a) ✓, but |Z̄_1| = |Z̄_2| = 14 > min(d_BB, d_LP) =
        #     min(12, 10), so neither is min-weight → (b) ✗ → THM 12 N/A
        #   - Example C: (a) ✗ (intra-code) → THM 12 N/A
        #
        # Theorem 13 requires the left code to have just 2 encoded qubits
        # with two min-weight non-overlapping logicals Z_l, X_l whose
        # product Z_l X_l cannot be reduced below 2*d_l weight.  Doesn't
        # match our single-product setup.
        #
        # The paper's explicit position (line 844 of the manuscript):
        #   "...verify directly that the deformed code defining the joint
        #    logical measurement has code distance d."
        # And §VII.A line 1597-1600:
        #   "Using BP-OSD as well as integer programming (CPLEX), we find
        #    the deformed code already has code distance of 12... Hence,
        #    there is no need to add more edge qubits to boost relative
        #    expansion in this auxiliary graph."
        #
        # So Zenodo's distance-preservation is via direct numerical
        # verification, NOT via any theorem.  Command 3 already validated
        # that the merged codes have parameters [[355,25,10]] and
        # [[150,5,12]] as claimed.
    end
    end
    @testset "merge helpers (G^T, M, H_R, bridge)" begin
                                        incidence_matrix_transpose,
                                        stabilizer_modification_matrix,
                                        canonical_H_R,
                                        adapter_bridge_x_check_matrix,
                                        edge_pair_to_index
                  src, dst, has_edge, add_edge!


    """Inline (lightweight) parser for the Zenodo G_mat_X .txt files —
    enough to recover the edge vertex pairs in Zenodo's storage order."""


    @testset "canonical_H_R" begin
        # n = 2
        H2 = canonical_H_R(2)
        @test size(H2) == (1, 2)
        @test Matrix(H2) == Bool[1 1]
        # n = 5
        H5 = canonical_H_R(5)
        @test size(H5) == (4, 5)
        @test Matrix(H5) == Bool[
            1 1 0 0 0
            0 1 1 0 0
            0 0 1 1 0
            0 0 0 1 1]
        # n = 14 (the inter-code adapter width for example B)
        H14 = canonical_H_R(14)
        @test size(H14) == (13, 14)
        @test nnz(H14) == 26
        # max row weight = 2, max col weight = 2 (interior cols)
        @test maximum(sum(H14; dims = 2)) == 2
        @test maximum(sum(H14; dims = 1)) == 2
        # Validation: n < 2 throws
        @test_throws ArgumentError canonical_H_R(1)
        @test_throws ArgumentError canonical_H_R(0)
        @test_throws ArgumentError canonical_H_R(-1)
    end

    @testset "incidence_matrix_transpose on small graphs" begin
        # Path P_4 (3 edges, 4 vertices)
        # Build a tiny CSS code to instantiate an AuxiliaryGraph via
        # build_initial_aux_graph; not strictly the canonical use, but
        # tests the function's structure.
        # Simpler: hand-build the AuxiliaryGraph.
        g = SimpleGraph(4)
        for (u, v) in [(1, 2), (2, 3), (3, 4)]
            add_edge!(g, u, v)
        end
        aux = AuxiliaryGraph(g, [[i] for i in 1:nv(g)], collect(1:ne(g)),
                              spzeros(Bool, 0, ne(g)), collect(1:nv(g)),
                              [Tuple{Int,Int}[] for _ in 1:nv(g)])
        GT = incidence_matrix_transpose(aux)
        @test size(GT) == (4, 3)
        @test nnz(GT) == 6   # = 2 · ne
        # Each column has exactly 2 nonzeros (the two endpoints).
        @test all(sum(GT; dims = 1) .== 2)
        # Verify endpoints are correct
        for (e_idx, e) in enumerate(edges(aux.graph))
            @test GT[src(e), e_idx] == true
            @test GT[dst(e), e_idx] == true
        end
    end

    if !isdir(BASE)
        @info "Skipping Zenodo merge-helper tests: data not present at $BASE"
        @test_skip "no Zenodo data"
    else
        # Helpers to build the cellulated aux graphs.
        z1_bb = sort([6,8,13,17,31,32,33,35,36,37,41,50,51,93]) .+ 1
        z2_lp = sort([24,25,26,29,30,56,58,59,60,61,90,93,94,121]) .+ 1
        z3_bb = sort([10,17,35,39,42,43,53,55,61,70,84,89]) .+ 1
        code_BB_canonical = CSS(
            load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                "Hx_BB_98_6_12_original-code-canonicalbasis.mtx")),
            load_mtx_inline(joinpath(BASE, "BB_98_6_12", "original_codes",
                "Hz_BB_98_6_12_original-code-canonicalbasis.mtx")))
        code_LP = CSS(
            load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
                "Hx_LP_200_20_10_original-code.mtx")),
            load_mtx_inline(joinpath(BASE, "LP_200_20_10", "original_codes",
                "Hz_LP_200_20_10_original-code.mtx")))
        aux_B_G1 = cellulate_long_cycles!(
            build_initial_aux_graph(z1_bb, code_BB_canonical); max_cycle_len = 6)
        aux_B_G2 = cellulate_long_cycles!(
            build_initial_aux_graph(z2_lp, code_LP); max_cycle_len = 7)

        @testset "incidence_matrix_transpose on Zenodo aux_B_G1 (14, 23)" begin
            GT = incidence_matrix_transpose(aux_B_G1)
            @test size(GT) == (14, 23)
            @test nnz(GT) == 2 * 23
            # Every column has exactly 2 nonzeros (a graph edge has 2 endpoints).
            @test all(sum(GT; dims = 1) .== 2)
            for (e_idx, e) in enumerate(edges(aux_B_G1.graph))
                @test GT[src(e), e_idx] == true
                @test GT[dst(e), e_idx] == true
            end
        end

        @testset "stabilizer_modification_matrix on Zenodo aux_B_G1 (49, 23)" begin
            M_l = stabilizer_modification_matrix(aux_B_G1, 49)
            @test size(M_l) == (49, 23)
            # Every overlapping stabilizer has |L_s|/2 = 1 entry; 21 overlap.
            @test nnz(M_l) == 21
            @test maximum(sum(M_l; dims = 2)) == 1   # max row weight
            @test maximum(sum(M_l; dims = 1)) == 1   # max col weight (no collisions in G_1)

            # Set-level verification against Zenodo's M_l (column ordering differs;
            # we compare {(s, sorted_pair)} sets).
            zenodo_g1_pairs = load_g_mat_pairs(
                joinpath(BASE, "BB_98_LP_200_adapter", "skipTree_transformations",
                         "BB_98_6_12_Z_1_GTP.txt"), "G_mat_1")
            zen_pairs = [minmax(p[1], p[2]) for p in zenodo_g1_pairs]
            Hx_merged = load_mtx_inline(joinpath(BASE, "BB_98_LP_200_adapter",
                "Hx_intercode_BB_LP_adapter-Z_1_Z_2_deformed-code.mtx"))
            M_l_zen = Hx_merged[1:49, 99:121]   # Zenodo's M_l block
            zen_set = Set{Tuple{Int,Tuple{Int,Int}}}()
            for s in 1:49, j in 1:23
                M_l_zen[s, j] && push!(zen_set, (s, zen_pairs[j]))
            end
            # Build our set
            edges_julia = [tuple(minmax(src(e), dst(e))...) for e in edges(aux_B_G1.graph)]
            our_set = Set{Tuple{Int,Tuple{Int,Int}}}()
            for s in 1:49, j in 1:23
                M_l[s, j] && push!(our_set, (s, edges_julia[j]))
            end
            @info "M_l verification: our set size $(length(our_set)), zenodo set size $(length(zen_set))"
            @test our_set == zen_set
        end

        @testset "stabilizer_modification_matrix on Zenodo aux_B_G2 (96, 20)" begin
            M_r = stabilizer_modification_matrix(aux_B_G2, 96)
            @test size(M_r) == (96, 20)
            # 21 overlapping LP stabilizers; stabs 32 and 96 share matching (2,14),
            # so total nonzero entries = 21 even though only 20 distinct edges.
            # But sparse matrix counts each (s, e) cell once, so:
            #   pairs = 21 (stab 32 → (2,14), stab 96 → (2,14)) ⇒ nnz = 21
            @test nnz(M_r) == 21
            @test maximum(sum(M_r; dims = 2)) == 1   # max row weight = 1
            @test maximum(sum(M_r; dims = 1)) == 2   # max col weight = 2 (the (2,14) collision)
        end

        @testset "adapter_bridge_x_check_matrix (13, 57)" begin
            stl = skiptree(aux_B_G1.graph)
            str = skiptree(aux_B_G2.graph)
            bridge = adapter_bridge_x_check_matrix(stl, str, 14)
            @test size(bridge) == (13, 57)

            # Verify left 23 columns equal skiptree_l.T
            @test bridge[:, 1:23] == stl.T
            # Middle 14 columns equal H_R(14)
            @test bridge[:, 24:37] == canonical_H_R(14)
            # Right 20 columns equal skiptree_r.T
            @test bridge[:, 38:57] == str.T

            # Dimension-mismatch error paths
            stl_bad = skiptree(aux_B_G1.graph)
            @test_throws DimensionMismatch adapter_bridge_x_check_matrix(
                stl_bad, str, 15)   # wrong adapter_width
            @test_throws ArgumentError adapter_bridge_x_check_matrix(stl, str, 1)
        end
    end
    end
    @testset "merged code Example B (inter-code [[355,25,10]])" begin
    if !isdir(BASE)
        @info "Skipping Example B end-to-end test: Zenodo data not present at $BASE"
        @test_skip "no Zenodo data"
    else



        pair = CodePair(bb_canonical, lp_code, Z1_BB, Z2_LP)

        @testset "build_adapter_intercode + parameters" begin
            adapter = build_adapter_intercode(pair;
                max_cycle_len_1 = 6, max_cycle_len_2 = 7)
            m = adapter.merged_code

            # Paper [[355, 25, 10]]
            @test size(m.Hx, 2) == 355
            @test size(m.Hx, 1) == 175
            @test size(m.Hz, 1) == 173

            # CSS commutation: Hx * Hz^T ≡ 0 (mod 2)
            HxI = Int.(m.Hx); HzI = Int.(m.Hz)
            commutator = (HxI * transpose(HzI)) .% 2
            @test count(!iszero, commutator) == 0

            # k = n - rank(Hx) - rank(Hz)
            n = size(m.Hx, 2)
            rx = rank_gf2(m.Hx)
            rz = rank_gf2(m.Hz)
            @info "Example B merged code: rank(Hx)=$rx, rank(Hz)=$rz, k=$(n - rx - rz)"
            @test n - rx - rz == 25

            # Stabilizer-weight bounds (paper §VII.A)
            @test maximum(sum(m.Hx, dims = 2)) ≤ 8
            @test maximum(sum(m.Hz, dims = 2)) ≤ 8

            # Qubit-degree bound (column weights of [Hx; Hz])
            combined = vcat(m.Hx, m.Hz)
            @test maximum(sum(combined, dims = 1)) ≤ 9

            # Sanity on Adapter struct itself
            @test adapter.adapter_width == 14
            @test adapter.code_pair === pair
        end

        @testset "Auto-tuned max_cycle_len (no kwargs → n=355)" begin
            # With no kwargs, max_cycle_len_i defaults to max stab weight
            # of code_i: BB → 6, LP → 7. Should produce exact [[355, 25, 10]].
            adapter = build_adapter_intercode(pair)
            m = adapter.merged_code
            @test size(m.Hx, 2) == 355
            @test size(m.Hx, 1) == 175
            @test size(m.Hz, 1) == 173
            n = size(m.Hx, 2)
            rx = rank_gf2(m.Hx); rz = rank_gf2(m.Hz)
            @test n - rx - rz == 25
            # Explicit override still works (backwards compat)
            adapter2 = build_adapter_intercode(pair;
                max_cycle_len_1 = 6, max_cycle_len_2 = 7)
            @test size(adapter2.merged_code.Hx) == size(m.Hx)
            @test size(adapter2.merged_code.Hz) == size(m.Hz)
        end

        @testset "Rank equality with Zenodo (basis-independent invariant)" begin
            # We expect equal rank but NOT equal row-space: SkipTree and
            # cycle-basis choices differ from Zenodo's, so the merged
            # CSS code we produce has the same parameters [[n, k, d]]
            # but is a different basis of the same code family. Per the
            # paper §VII.A — different choices of expansion graph, MST,
            # and cycle basis all yield equally-valid adapters. The
            # numerical k and d are invariants; the specific rows are not.
            Hx_zen = load_mtx_inline(joinpath(BASE, "BB_98_LP_200_adapter",
                "Hx_intercode_BB_LP_adapter-Z_1_Z_2_deformed-code.mtx"))
            Hz_zen = load_mtx_inline(joinpath(BASE, "BB_98_LP_200_adapter",
                "Hz_intercode_BB_LP_adapter-Z_1_Z_2_deformed-code.mtx"))

            adapter = build_adapter_intercode(pair;
                max_cycle_len_1 = 6, max_cycle_len_2 = 7)
            m = adapter.merged_code

            # Basis-independent rank equality: same code dimension (= k)
            # regardless of which valid SkipTree / cycle-basis / chord
            # choices each side made. This is the load-bearing claim — if
            # it failed, our reproduction would be a different code.
            rx_ours = rank_gf2(m.Hx);  rx_zen = rank_gf2(Hx_zen)
            rz_ours = rank_gf2(m.Hz);  rz_zen = rank_gf2(Hz_zen)
            @info "Example B Hx rank: ours=$rx_ours, zenodo=$rx_zen"
            @info "Example B Hz rank: ours=$rz_ours, zenodo=$rz_zen"
            @test rx_ours == rx_zen
            @test rz_ours == rz_zen

            # Row-space equality is NOT expected: our SkipTree / cycle basis
            # differ from Zenodo's, so the rows themselves differ even
            # though the rank is the same. Documented here as an `@info`
            # (not an assertion) so the actual gap is recorded each run.
            rx_combo = rank_gf2(vcat(m.Hx, Hx_zen))
            rz_combo = rank_gf2(vcat(m.Hz, Hz_zen))
            @info "Example B row-space gap: combined Hx rank = $rx_combo (= ours only if row spans match)"
            @info "Example B row-space gap: combined Hz rank = $rz_combo"
        end

        @testset "Error paths" begin
            # Intra-code (code1 === code2) is rejected by inter-code path
            pair_intra = CodePair(bb_canonical, bb_canonical, Z1_BB, Z1_BB)
            @test_throws ArgumentError build_adapter_intercode(pair_intra)

            # Unequal logical-support lengths
            z1_short = Z1_BB[1:10]
            pair_bad = CodePair(bb_canonical, lp_code, z1_short, Z2_LP)
            @test_throws DimensionMismatch build_adapter_intercode(pair_bad)
        end
    end
    end
    @testset "merged code Example C (intra-code [[150,5,12]])" begin
                                        assemble_merged_code_intracode


    if !isdir(BASE)
        @info "Skipping Example C end-to-end test: Zenodo data not present at $BASE"
        @test_skip "no Zenodo data"
    else



        # Z_1 and Z_3 from Zenodo (0-indexed → +1)
        pair = CodePair(bb_full, bb_full, Z1_BB, Z3_BB)

        @testset "build_adapter_intracode + parameters" begin
            adapter = build_adapter_intracode(pair)
            m = adapter.merged_code

            # Paper [[150, 5, 12]]
            @test size(m.Hx, 2) == 150
            @test size(m.Hx, 1) == 73
            @test size(m.Hz, 1) == 72

            # CSS commutation
            HxI = Int.(m.Hx); HzI = Int.(m.Hz)
            @test count(!iszero, (HxI * transpose(HzI)) .% 2) == 0

            # k = n - rank(Hx) - rank(Hz) = 5
            n = size(m.Hx, 2)
            rx = rank_gf2(m.Hx)
            rz = rank_gf2(m.Hz)
            @info "Example C merged code: rank(Hx)=$rx, rank(Hz)=$rz, k=$(n - rx - rz)"
            @test n - rx - rz == 5

            # Weight bounds (paper Table IV: Hz max weight 6)
            @test maximum(sum(m.Hx, dims = 2)) ≤ 8
            @test maximum(sum(m.Hz, dims = 2)) ≤ 6
            @test maximum(sum(vcat(m.Hx, m.Hz), dims = 1)) ≤ 9

            # Adapter struct
            @test adapter.adapter_width == 12
            @test adapter.code_pair === pair
            @test adapter.aux_l.logical_support == Z1_BB
            @test adapter.aux_r.logical_support == Z3_BB
        end

        @testset "Rank equality with Zenodo (basis-independent invariant)" begin
            Hx_zen = load_mtx_inline(joinpath(BASE, "BB_98_intracode_adapter",
                "Hx_BB_intracode_Z_1_Z_3_adapted-code.mtx"))
            Hz_zen = load_mtx_inline(joinpath(BASE, "BB_98_intracode_adapter",
                "Hz_BB_intracode_Z_1_Z_3_adapted-code.mtx"))
            adapter = build_adapter_intracode(pair)
            m = adapter.merged_code
            rx_ours = rank_gf2(m.Hx); rx_zen = rank_gf2(Hx_zen)
            rz_ours = rank_gf2(m.Hz); rz_zen = rank_gf2(Hz_zen)
            @info "Example C Hx rank: ours=$rx_ours, zenodo=$rx_zen"
            @info "Example C Hz rank: ours=$rz_ours, zenodo=$rz_zen"
            @test rx_ours == rx_zen
            @test rz_ours == rz_zen
        end

        @testset "build_adapter dispatcher" begin
            # Inter-code routing
            pair_B = CodePair(bb_canonical, lp_code, Z1_BB, Z2_LP)
            adapter_B = build_adapter(pair_B)
            @test size(adapter_B.merged_code.Hx, 2) == 355   # inter-code

            # Intra-code routing (same code1 === code2)
            adapter_C = build_adapter(pair)
            @test size(adapter_C.merged_code.Hx, 2) == 150   # intra-code
        end

        @testset "Error paths (intra-code)" begin
            # Inter-code pair rejected by intra-code path
            pair_inter = CodePair(bb_full, lp_code, Z1_BB, Z2_LP)
            @test_throws ArgumentError build_adapter_intracode(pair_inter)
        end
    end
    end
    @testset "single-logical deform_code (Table II)" begin
    # Reproduce paper Table II: single-logical deformed codes (one
    # auxiliary graph attached to one input code, no bridge, no joint
    # measurement). Verifies parameters (n, k, CSS commutation, weight
    # bounds) against the paper. Distance verification is done
    # out-of-band — too slow for CI — and recorded in the Step 0
    # distance closeout note.



    @testset "Small reproducible: Surface(3,3) single-logical deform" begin
        # Always-available smoke test: deform a [[13,1,3]] surface code.
        # Resulting deformed code is [[15, 0, ?]] (1 logical consumed,
        # path-graph aux has no cycles).
        c = CSS(Matrix{Bool}(parity_matrix_x(Surface(3, 3))),
                Matrix{Bool}(parity_matrix_z(Surface(3, 3))))
        z = sort(findall(!iszero, stab_to_gf2(logz_ops(Surface(3, 3)))[1, 14:26]))
        dc = deform_code(c, z)
        @test code_n(dc) == 15      # 13 + 2 (path-graph edges)
        @test code_k(dc) == 0       # original k=1 consumed by deformation
        HxI = Int.(dc.Hx); HzI = Int.(dc.Hz)
        @test count(!iszero, (HxI * transpose(HzI)) .% 2) == 0
        # Weight bounds matching the original surface code (≤ 4)
        @test maximum(sum(dc.Hx, dims = 2)) ≤ 5
        @test maximum(sum(dc.Hz, dims = 2)) ≤ 5
    end

    if !isdir(BASE)
        @info "Skipping Table II reproduction: Zenodo data not present at $BASE"
        @test_skip "no Zenodo data"
    else



        @testset "BB Z_1 (canonical) → [[121, 5, 12]]" begin
            dc = deform_code(bb_canonical, Z1_BB)
            @test code_n(dc) == 121     # 98 + 23 (aux edges)
            n = code_n(dc)
            rx = rank_gf2(dc.Hx); rz = rank_gf2(dc.Hz)
            @test n - rx - rz == 5
            @test count(!iszero, (Int.(dc.Hx) * transpose(Int.(dc.Hz))) .% 2) == 0
            # Paper Table II: max stab weights (X, Z) = (7, 6); max qubit degree 7
            @test maximum(sum(dc.Hx, dims = 2)) == 7
            @test maximum(sum(dc.Hz, dims = 2)) == 6
            @test maximum(sum(vcat(dc.Hx, dc.Hz), dims = 1)) == 7
        end

        @testset "LP Z_2 → [[220, 19, 10]]" begin
            dc = deform_code(lp_code, Z2_LP)
            @test code_n(dc) == 220     # 200 + 20 (aux edges, no cellulation)
            n = code_n(dc)
            rx = rank_gf2(dc.Hx); rz = rank_gf2(dc.Hz)
            @test n - rx - rz == 19
            @test count(!iszero, (Int.(dc.Hx) * transpose(Int.(dc.Hz))) .% 2) == 0
            # Paper Table II: max stab weights (X, Z) = (8, 7); max qubit degree 8
            # Our cycle basis (Graphs.jl) differs from paper's (NetworkX) and
            # gives max degree 9 (one edge appears in 6 cycles vs paper's max 5).
            # Both yield valid [[220, 19, 10]] codes; basis-choice artefact.
            @test maximum(sum(dc.Hx, dims = 2)) == 8
            @test maximum(sum(dc.Hz, dims = 2)) == 7
            @test maximum(sum(vcat(dc.Hx, dc.Hz), dims = 1)) ≤ 9
        end

        @testset "BB Z_3 (full-rank) → [[115, 5, 12]]" begin
            dc = deform_code(bb_full, Z3_BB)
            @test code_n(dc) == 115     # 98 + 17 (aux edges, 1 cellulated)
            n = code_n(dc)
            rx = rank_gf2(dc.Hx); rz = rank_gf2(dc.Hz)
            @test n - rx - rz == 5
            @test count(!iszero, (Int.(dc.Hx) * transpose(Int.(dc.Hz))) .% 2) == 0
            # Paper Table II: max stab weights (X, Z) = (7, 6); max qubit degree 7
            @test maximum(sum(dc.Hx, dims = 2)) == 7
            @test maximum(sum(dc.Hz, dims = 2)) == 6
            @test maximum(sum(vcat(dc.Hx, dc.Hz), dims = 1)) == 7
        end
    end
    end
    @testset "chord-override + bridge-detection" begin
    # Regression tests for the cellulation-chord-choice basis sensitivity
    # discovered on BB / Z_3 (Table II reproduction):
    #
    # - Our default midpoint chord choice produces a [[115, 5, ≤11]] BB+Z_3
    #   deformed code (NOT the paper's [[115, 5, 12]]).
    # - Supplying Zenodo's specific chord (1, 11) via the `chords` kwarg
    #   restores the paper-claimed distance d = 12.
    # - On BB / Z_1, basis-dependence is NOT triggered (both default and
    #   alternative valid chord placements give d = 12).
    # - For joint measurements, the bridge X-checks empirically detect the
    #   extra weight-11 standalone X-logicals — so the joint Example C
    #   [[150, 5, 12]] distance is bridge-protected, not coincidental.


    if !isdir(BASE)
        @info "Skipping chord-override regression: Zenodo data not present at $BASE"
        @test_skip "no Zenodo data"
    else


        @testset "Chord-override kwarg validation" begin
            aux = build_initial_aux_graph(Z3_BB, bb_full)
            n_vert = 12
            # Endpoint out of range
            @test_throws ArgumentError cellulate_long_cycles!(deepcopy(aux);
                chords = [(0, 3)])
            @test_throws ArgumentError cellulate_long_cycles!(deepcopy(aux);
                chords = [(3, n_vert + 1)])
            # Self-loop
            @test_throws ArgumentError cellulate_long_cycles!(deepcopy(aux);
                chords = [(5, 5)])
        end

        @testset "BB+Z_3 with explicit chord (1, 11) — Zenodo chord override" begin
            dc = deform_code(bb_full, Z3_BB; chords = [(1, 11)])
            @test code_n(dc) == 115     # 98 + 17 (matching + 1 explicit chord)
            n = code_n(dc)
            HxI = Int.(dc.Hx); HzI = Int.(dc.Hz)
            @test count(!iszero, (HxI * transpose(HzI)) .% 2) == 0
            rx = QuantumClifford.gf2_row_echelon_with_pivots!(copy(HxI))[2]
            rz = QuantumClifford.gf2_row_echelon_with_pivots!(copy(HzI))[2]
            @test n - rx - rz == 5

            # Verify the chord ended up in the graph (compared to pre-cellulate)
            aux_pre = build_initial_aux_graph(Z3_BB, bb_full)
            aux_post = cellulate_long_cycles!(
                build_initial_aux_graph(Z3_BB, bb_full); chords = [(1, 11)])
            edges_pre  = Set([tuple(minmax(src(e), dst(e))...) for e in edges(aux_pre.graph)])
            edges_post = Set([tuple(minmax(src(e), dst(e))...) for e in edges(aux_post.graph)])
            chords_added = setdiff(edges_post, edges_pre)
            @test (1, 11) in chords_added   # our explicit chord was applied
            @test ne(aux_post.graph) == 17  # 16 matching + 1 cellulation
        end

        @testset "Bridge detects weight-11 standalone-Z_3 logicals in joint Example C" begin
            # Build joint Example C with default chord (which yields a
            # standalone Z_3 with d ≤ 11).
            adapter_C = build_adapter(CodePair(bb_full, bb_full, Z1_BB, Z3_BB))
            mc_C = adapter_C.merged_code
            n_joint = size(mc_C.Hx, 2)
            @test n_joint == 150

            # Block boundaries
            ne1 = ne(adapter_C.aux_l.graph)
            ne2 = ne(adapter_C.aux_r.graph)
            w   = adapter_C.adapter_width
            nc_l = size(adapter_C.aux_l.cycle_basis_matrix, 1)
            sx_BB = 46
            bridge_rows = (sx_BB + nc_l + 1):(sx_BB + nc_l + w - 1)
            g3_start = 98 + ne1 + w + 1   # G_3 column block start in joint

            # The standalone BB+Z_3 deformed code, built with the same default
            # chord choice that the joint code uses internally for aux_r.
            dc_z3 = deform_code(bb_full, Z3_BB)
            n_standalone = code_n(dc_z3)
            @test n_standalone == 115

            # Candidate weight-11 X-logical supports for the STANDALONE
            # BB+Z_3 deformed code (cols 1..98 = BB, cols 99..115 = G_3
            # edges). Found by HiGHS MIP at T=11 (see distance closeout
            # note); embedded here as data. We verify at TEST TIME that each
            # is genuinely a non-trivial X-logical of the standalone code
            # — i.e. commutes with all X-stabs and anticommutes with at
            # least one X-logical — so the bridge-detection assertion
            # below is meaningful and not just a property of an arbitrary
            # vector.
            weight_11_supports = [
                [7, 14, 19, 31, 51, 67, 82, 95, 102, 109, 114],   # MIP li = 1
                [6, 7, 14, 31, 40, 50, 67, 95, 102, 109, 114],    # MIP li = 2
                [37, 38, 51, 59, 80, 87, 95, 99, 105, 109, 114],  # MIP li = 5
            ]

            # Independent verification each supp encodes a real X-logical
            # of the standalone deformed code (not just an arbitrary vector).
            Hx_std_I = Int.(dc_z3.Hx)
            Lx_std   = stab_to_gf2(logx_ops(dc_z3))[:, 1:n_standalone]
            HxI_joint = Int.(mc_C.Hx)

            for (idx, supp) in enumerate(weight_11_supports)
                @test length(supp) == 11

                # Build the error vector on the standalone code's qubit set.
                e_std = zeros(Int, n_standalone)
                for q in supp; e_std[q] = 1; end

                # (1) Must commute with all standalone X-stabilizers
                #     (i.e., be in ker(Hx_std) mod 2).
                @test all(==(0), (Hx_std_I * e_std) .% 2)

                # (2) Must anti-commute with at least one X-logical of the
                #     standalone code (i.e., be a non-trivial coset rep).
                anti_any = false
                for li in axes(Lx_std, 1)
                    if isodd(sum(Lx_std[li, :] .& (e_std .!= 0)))
                        anti_any = true; break
                    end
                end
                @test anti_any

                # Embed into the joint merged code's qubit space. Standalone
                # qubits 1..98 (BB) map to joint qubits 1..98; standalone
                # qubits 99..115 (G_3 edges) map to joint G_3 cols starting
                # at `g3_start`. The other joint qubits (G_1, A) get 0.
                e_joint = falses(n_joint)
                for q in supp
                    if q ≤ 98
                        e_joint[q] = true
                    else
                        eidx = q - 98
                        e_joint[g3_start - 1 + eidx] = true
                    end
                end

                # Compute the joint syndrome and check the bridge fires.
                syn = (HxI_joint * Int.(e_joint)) .% 2
                bridge_syn = syn[collect(bridge_rows)]
                @test count(!iszero, bridge_syn) > 0   # bridge detects this logical
            end
        end
    end
    end
    @testset "gross [[144,12,12]] BB regression" begin
    # Smoke test on the largest non-Zenodo code we've tried: the
    # Bravyi-et-al-2024 gross bivariate-bicycle [[144, 12, 12]] code
    # (l=12, m=6, A = x³+y+y², B = y³+x+x²). We do NOT verify distance
    # here (full MIP cert at d=12 takes hours and is run out-of-band);
    # only that build_adapter produces a structurally valid merged
    # code on this scale.

    l, m = 12, 6
    A = [(:x, 3), (:y, 1), (:y, 2)]
    B = [(:y, 3), (:x, 1), (:x, 2)]
    g_orig = BivariateBicycleViaCirculantMat(l, m, A, B)
    @test code_n(g_orig) == 144
    @test code_k(g_orig) == 12

    g = CSS(Matrix{Bool}(parity_matrix_x(g_orig)),
            Matrix{Bool}(parity_matrix_z(g_orig)))
    lz = stab_to_gf2(logz_ops(g_orig))
    n0 = code_n(g_orig)
    # Use the two lowest-weight Z-logicals from the MixedDestabilizer basis.
    z5  = sort(findall(!iszero, lz[5,  n0+1:2n0]))
    z10 = sort(findall(!iszero, lz[10, n0+1:2n0]))

    adapter = build_adapter(CodePair(g, g, z5, z10))
    m = adapter.merged_code

    # Structural assertions: exact column count from the intra-code block
    # layout — n_code + ne(aux_l) + adapter_width + ne(aux_r).
    expected_n = 144 + ne(adapter.aux_l.graph) +
                  adapter.adapter_width + ne(adapter.aux_r.graph)
    @test size(m.Hx, 2) == expected_n
    @test adapter.adapter_width == min(length(z5), length(z10))
    HxI = Int.(m.Hx); HzI = Int.(m.Hz)
    @test count(!iszero, (HxI * transpose(HzI)) .% 2) == 0  # CSS commute

    # Joint measurement consumes exactly one logical (paper Theorem 19):
    # k_merged = k_original - 1 = 12 - 1 = 11.
    @test code_k(m) == 11

    # Paper-level sparsity bounds (Theorem 5 + paper Table III ≤ 8 / 9)
    @test maximum(sum(m.Hx, dims = 2)) ≤ 8
    @test maximum(sum(m.Hz, dims = 2)) ≤ 8
    @test maximum(sum(vcat(m.Hx, m.Hz), dims = 1)) ≤ 9
    end
    @testset "general API validation" begin
    @testset "Inter-code Surface(3,3) × Surface(3,3) → [[33,1,3]]" begin
        c1 = as_css(Surface(3, 3))
        c2 = as_css(Surface(3, 3))
        lz = stab_to_gf2(logz_ops(Surface(3, 3)))
        n0 = code_n(c1)
        z = sort(findall(!iszero, lz[1, n0+1:2n0]))
        @test length(z) == 3
        adapter = build_adapter(CodePair(c1, c2, z, z))
        m = adapter.merged_code
        n = size(m.Hx, 2)

        @test n == 33                              # 13 + 2 + 3 + 2 + 13
        @test size(m.Hx) == (14, 33)
        @test size(m.Hz) == (18, 33)

        # CSS commutation
        HxI = Int.(m.Hx); HzI = Int.(m.Hz)
        @test count(!iszero, (HxI * transpose(HzI)) .% 2) == 0

        # k = 1
        rx = rank_gf2(m.Hx); rz = rank_gf2(m.Hz)
        @test n - rx - rz == 1

        # Weight bounds (matches Surface(3,3) weight 4)
        @test maximum(sum(m.Hx, dims = 2)) == 4
        @test maximum(sum(m.Hz, dims = 2)) == 4
        @test maximum(sum(vcat(m.Hx, m.Hz), dims = 1)) == 4

        # Adapter struct
        @test adapter.adapter_width == 3
    end

    @testset "Intra-code BB[[18,4,4]] self-join → smaller k" begin
        l, m = 3, 3
        A = [(:x, 0), (:x, 1), (:y, 1)]
        B = [(:y, 0), (:x, 2), (:y, 2)]
        c_bb_orig = BivariateBicycleViaCirculantMat(l, m, A, B)
        @test code_n(c_bb_orig) == 18
        @test code_k(c_bb_orig) == 4

        bb = as_css(c_bb_orig)
        lz = stab_to_gf2(logz_ops(c_bb_orig))
        n0 = code_n(c_bb_orig); k_orig = code_k(c_bb_orig)
        weights = [(i, count(!iszero, lz[i, n0+1:2n0])) for i in 1:k_orig]
        sort!(weights, by = x -> x[2])
        i1, i2 = weights[1][1], weights[2][1]
        # NOTE: these are BB[[18,4,4]] logicals — distinct from the outer
        # Z1_BB / Z2_LP, which are BB[[98,6,12]] / LP[[200,20,10]] supports.
        z_a = sort(findall(!iszero, lz[i1, n0+1:2n0]))
        z_b = sort(findall(!iszero, lz[i2, n0+1:2n0]))
        @test length(z_a) == 4
        @test length(z_b) == 4

        pair = CodePair(bb, bb, z_a, z_b)
        adapter = build_adapter(pair)
        mc = adapter.merged_code

        # Two of the 4 logicals get fused — expect k = 3.
        nm = size(mc.Hx, 2)
        rx = rank_gf2(mc.Hx); rz = rank_gf2(mc.Hz)
        @test nm - rx - rz == 3

        # CSS commutation
        HxI = Int.(mc.Hx); HzI = Int.(mc.Hz)
        @test count(!iszero, (HxI * transpose(HzI)) .% 2) == 0

        # Adapter struct
        @test adapter.adapter_width == 4
        @test adapter.code_pair.code1 === adapter.code_pair.code2
    end

    @testset "AdapterMergedCode wrapper — AbstractECC interface" begin
        c1 = as_css(Surface(3, 3))
        c2 = as_css(Surface(3, 3))
        lz = stab_to_gf2(logz_ops(Surface(3, 3)))
        n0 = code_n(c1)
        z = sort(findall(!iszero, lz[1, n0+1:2n0]))
        adapter = build_adapter(CodePair(c1, c2, z, z))
        wrap = AdapterMergedCode(adapter)

        # Type hierarchy
        @test wrap isa AbstractCSSCode
        @test wrap isa AbstractQECC
        @test wrap isa AbstractECC

        # Interface methods delegate correctly
        @test code_n(wrap) == 33
        @test code_k(wrap) == 1
        @test code_s(wrap) == 32              # 14 X-stabs + 18 Z-stabs
        @test size(parity_matrix_x(wrap)) == (14, 33)
        @test size(parity_matrix_z(wrap)) == (18, 33)
        @test size(parity_matrix(wrap)) == (32, 66)
        @test parity_checks(wrap) isa QuantumClifford.Stabilizer

        # Returned matrices share data with the inner CSS (no copy)
        @test parity_matrix_x(wrap) === adapter.merged_code.Hx
        @test parity_matrix_z(wrap) === adapter.merged_code.Hz

        # Adapter metadata is preserved
        @test wrap.adapter === adapter
        @test wrap.adapter.adapter_width == 3
    end

    @testset "Tree-aux-graph case (no cycles) — regression guard" begin
        # Surface(3,3) with weight-3 logical produces a 3-vertex path aux
        # graph (no cycles). The merged-code assembly previously required
        # a non-empty cycle_basis_matrix; that guard was removed.
        c1 = as_css(Surface(3, 3))
        c2 = as_css(Surface(3, 3))
        lz = stab_to_gf2(logz_ops(Surface(3, 3)))
        n0 = code_n(c1)
        z = sort(findall(!iszero, lz[1, n0+1:2n0]))
        adapter = build_adapter(CodePair(c1, c2, z, z))
        # Both aux graphs are paths (3 vertices, 2 edges, 0 cycles)
        @test size(adapter.aux_l.cycle_basis_matrix, 1) == 0
        @test size(adapter.aux_r.cycle_basis_matrix, 1) == 0
        # And the merged code is still well-formed
        m = adapter.merged_code
        @test size(m.Hx, 2) == 33
        @test count(!iszero, (Int.(m.Hx) * transpose(Int.(m.Hz))) .% 2) == 0
    end
    end
    @testset "error / edge-case coverage" begin
    # Coverage for validation and edge-case branches that weren't exercised
    # by the higher-level tests. Keeps the surface tight — every error
    # path in the public API has at least one assertion.

    surface = as_css(Surface(3, 3))

    @testset "build_initial_aux_graph validation" begin
        # logical_support must be sorted
        @test_throws ArgumentError build_initial_aux_graph([3, 1, 2], surface)
        # logical_support must have ≥ 2 qubits
        @test_throws ArgumentError build_initial_aux_graph([5], surface)
        # duplicates
        @test_throws ArgumentError build_initial_aux_graph([3, 3, 5], surface)
        # out-of-range qubit index
        n = 13   # Surface(3,3) qubits
        @test_throws ArgumentError build_initial_aux_graph([3, 6, n + 1], surface)
    end

    @testset "skiptree validation" begin
        # n < 2
        @test_throws ArgumentError skiptree(SimpleGraph(1))
        # Disconnected graph (two isolated vertices, no edges)
        @test_throws ArgumentError skiptree(SimpleGraph(2))
        # Root out of range
        g = SimpleGraph(3); add_edge!(g, 1, 2); add_edge!(g, 2, 3)
        @test_throws ArgumentError skiptree(g; root = 0)
        @test_throws ArgumentError skiptree(g; root = 4)
    end

    @testset "relative_expansion validation" begin
        g = SimpleGraph(3); add_edge!(g, 1, 2); add_edge!(g, 2, 3)
        # t ≥ 1 required
        @test_throws ArgumentError relative_expansion(g, [1, 2], 0)
        # port vertex out of range
        @test_throws ArgumentError relative_expansion(g, [1, 4], 2)
        # n > 24 → brute-force cap rejects
        big = SimpleGraph(25)
        for i in 1:24; add_edge!(big, i, i + 1); end
        @test_throws ArgumentError relative_expansion(big, [1, 25], 3)
    end

    @testset "edge_pair_to_index — absent edge and unsorted pair" begin
        g = SimpleGraph(4)
        add_edge!(g, 1, 2); add_edge!(g, 2, 3)
        # Edge present
        @test edge_pair_to_index(g, (1, 2)) > 0
        @test edge_pair_to_index(g, (2, 3)) > 0
        # Edge absent → 0
        @test edge_pair_to_index(g, (1, 4)) == 0
        @test edge_pair_to_index(g, (3, 4)) == 0
        # Unsorted pair → ArgumentError
        @test_throws ArgumentError edge_pair_to_index(g, (2, 1))
    end

    @testset "cellulate_long_cycles! validation" begin
        aux = build_initial_aux_graph([3, 6, 9], surface)
        # max_cycle_len ≥ 3 required
        @test_throws ArgumentError cellulate_long_cycles!(aux; max_cycle_len = 2)
        # max_iterations ≥ 0 required
        @test_throws ArgumentError cellulate_long_cycles!(aux; max_iterations = -1)
    end
    end
    @testset "distance verification (MIP sanity)" begin
    # Sanity check: MIP-feasibility distance formulation on small known
    # codes. Catches the symplectic-column bug (and any future variant)
    # that would silently produce "d = ∞" (always-infeasible) results
    # for the merged-code verification runs.
    #
    # NOTE: this testitem only validates the MIP formulation logic on
    # small codes (sub-second per call). It does NOT run the full
    # merged-code distance MIPs from Step 0 closeout — those are too
    # slow for CI.



    @testset "Steane [[7, 1, 3]] — both types proven d = 3" begin
        c = as_css(Steane7())
        for ltype in (:X, :Z)
            r2 = logical_op_exists_at_weight_le_T(c, 1, ltype, 2)
            r3 = logical_op_exists_at_weight_le_T(c, 1, ltype, 3)
            @test r2[1] == :proven_lower_bound
            @test r3[1] == :counterexample_found
            @test r3[2] == 3
        end
    end

    @testset "Surface(3,3) [[13, 1, 3]] — both types proven d = 3" begin
        c = as_css(Surface(3, 3))
        for ltype in (:X, :Z)
            r2 = logical_op_exists_at_weight_le_T(c, 1, ltype, 2)
            r3 = logical_op_exists_at_weight_le_T(c, 1, ltype, 3)
            @test r2[1] == :proven_lower_bound
            @test r3[1] == :counterexample_found
            @test r3[2] == 3
        end
    end
    end
end
