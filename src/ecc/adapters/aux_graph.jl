"""
Initial auxiliary graph G_0 from paper §II.C [swaroop2026universal](@cite).
One vertex per qubit in `logical_support`, one edge per pair of consecutive
vertices in `sort(supp(s) ∩ logical_support)` for each X-stabilizer that
overlaps the logical. Warns and adds repair edges if the result is disconnected.

Satisfies Theorem 5 desiderata 0–2 (connectivity, bounded degree, short
matchings). Desideratum 3 (sparse cycle basis) requires a follow-up call
to [`cellulate_long_cycles!`](@ref). Desideratum 4 (relative expansion)
is verified separately via [`verify_desideratum_4`](@ref).
"""
function build_initial_aux_graph(
    logical_support::Vector{Int},
    code::CSS,
)::AuxiliaryGraph
    n_support = length(logical_support)
    n_support ≥ 2 || throw(ArgumentError(
        "logical_support must have ≥ 2 qubits (got $n_support)"))

    issorted(logical_support) ||
        throw(ArgumentError("logical_support must be sorted"))
    allunique(logical_support) ||
        throw(ArgumentError("logical_support must have distinct qubit indices"))
    n_qubits = size(code.Hx, 2)
    for q in logical_support
        1 ≤ q ≤ n_qubits || throw(ArgumentError(
            "logical_support qubit $q is out of range 1..$n_qubits"))
    end

    qubit_to_vertex = Dict{Int,Int}(q => i for (i, q) in enumerate(logical_support))
    support_set = Set(logical_support)

    # Store stabilizer matchings as vertex pairs, not edge indices:
    # cellulation later adds chords and that shifts the edges() iteration order.
    g = SimpleGraph(n_support)
    n_stabs = size(code.Hx, 1)
    stabilizer_matchings = [Tuple{Int,Int}[] for _ in 1:n_stabs]

    @inbounds for s in 1:n_stabs
        supp_s = Int[]
        for q in 1:n_qubits
            if code.Hx[s, q] && q in support_set
                push!(supp_s, q)
            end
        end
        length(supp_s) < 2 && continue
        isodd(length(supp_s)) && error(
            "stabilizer $s overlaps logical on $(length(supp_s)) qubits (odd); " *
            "stabilizer must commute with the target Z̄ logical")
        vs = sort!([qubit_to_vertex[q] for q in supp_s])
        for k in 1:2:length(vs)
            v_a, v_b = vs[k], vs[k + 1]
            has_edge(g, v_a, v_b) || add_edge!(g, v_a, v_b)
            push!(stabilizer_matchings[s], (v_a, v_b))
        end
    end

    if !is_connected(g)
        @warn "auxiliary graph G_0 is disconnected after matching step; " *
              "adding repair edges (greedy nearest cross-component pairing)" *
              " — produced graph satisfies desideratum 0 but is not optimised."
        comps = connected_components(g)
        while length(comps) > 1
            c1 = comps[1]; c2 = comps[2]
            v_a, v_b = minmax(minimum(c1), minimum(c2))
            has_edge(g, v_a, v_b) || add_edge!(g, v_a, v_b)
            comps = connected_components(g)
        end
    end

    n_edges = ne(g)
    port_function = [[i] for i in 1:n_support]
    edge_to_qubit = collect(1:n_edges)   # placeholder; merge step rewrites
    cycle_basis_matrix = spzeros(Bool, 0, n_edges)

    AuxiliaryGraph(g, port_function, edge_to_qubit, cycle_basis_matrix,
                   copy(logical_support), stabilizer_matchings)
end

"""
The `F` matrix from Figure 7 [swaroop2026universal](@cite): `F[v, q] = 1`
iff vertex `v` is in the port function of qubit `q`. Shape
`(nv(aux.graph), n_qubits)`.
"""
function port_function_as_matrix(aux::AuxiliaryGraph,
                                  n_qubits::Int)::SparseMatrixCSC{Bool, Int}
    n_vert = nv(aux.graph)
    I = Int[]; J = Int[]
    for (q_idx, q) in enumerate(aux.logical_support)
        1 ≤ q ≤ n_qubits ||
            throw(ArgumentError("logical_support qubit $q out of range 1..$n_qubits"))
        for v in aux.port_function[q_idx]
            1 ≤ v ≤ n_vert ||
                throw(ArgumentError("vertex $v out of range 1..$n_vert"))
            push!(I, v); push!(J, q)
        end
    end
    sparse(I, J, trues(length(I)), n_vert, n_qubits)
end

"""
Position of the sorted vertex pair `(u, v)` in `edges(graph)` iteration order,
or `0` if it isn't an edge. This iteration order is the column convention
for the merged H_X / H_Z edge blocks.

O(`ne(graph)`). Hot loops should precompute `Dict{Tuple{Int,Int}, Int}`.
"""
function edge_pair_to_index(graph::SimpleGraph{Int},
                              pair::Tuple{Int,Int})::Int
    pair[1] ≤ pair[2] ||
        throw(ArgumentError("edge_pair_to_index: pair must be sorted (got $pair)"))
    has_edge(graph, pair[1], pair[2]) || return 0
    for (i, e) in enumerate(edges(graph))
        u, v = minmax(src(e), dst(e))
        if (u, v) == pair
            return i
        end
    end
    0  # unreachable if has_edge passed and graph is consistent
end

"""
Cellulate `aux.graph` so every basis cycle has length ≤ `max_cycle_len`
(paper desideratum 3a). Repeatedly takes the longest basis cycle and adds
the chord `(cycle[1], cycle[n÷2 + 1])`. If `chords` is supplied, those are
added in order first, then the midpoint algorithm runs.

Mutates `aux.graph` and `aux.edge_to_qubit`, rebuilds the cycle basis
matrix, and returns a new `AuxiliaryGraph` sharing the mutated graph.

# The `chords` kwarg — basis choice affects distance

`Graphs.cycle_basis` and NetworkX make different DFS choices, so the
midpoint cellulation differs from Zenodo's. For paper Table II's BB / Z_3
deformed code, the paper's chord `(1, 11)` gives `[[115, 5, 12]]` while
our default `(4, 9)` gives `[[115, 5, ≤11]]`. Pass `chords=[(1, 11)]` for
byte-identical paper reproduction. Joint measurements via
[`build_adapter`](@ref) are unaffected — the bridge X-checks absorb the
extra weight-`(d-1)` logical that the chord choice introduces.
"""
function cellulate_long_cycles!(
    aux::AuxiliaryGraph;
    max_cycle_len::Int = 6,
    max_iterations::Int = 100,
    chords::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
)::AuxiliaryGraph
    max_cycle_len ≥ 3 ||
        throw(ArgumentError("max_cycle_len must be ≥ 3 (got $max_cycle_len)"))
    max_iterations ≥ 0 ||
        throw(ArgumentError("max_iterations must be ≥ 0 (got $max_iterations)"))

    g = aux.graph
    n_vert = nv(g)

    if chords !== nothing
        for (i, (u, v)) in enumerate(chords)
            (1 ≤ u ≤ n_vert && 1 ≤ v ≤ n_vert) || throw(ArgumentError(
                "cellulate_long_cycles!: chord #$i ($u, $v) has endpoint out of range 1:$n_vert"))
            u != v || throw(ArgumentError(
                "cellulate_long_cycles!: chord #$i ($u, $v) is a self-loop"))
            !has_edge(g, u, v) || throw(ArgumentError(
                "cellulate_long_cycles!: chord #$i ($u, $v) already an edge"))
            add_edge!(g, u, v)
            push!(aux.edge_to_qubit, length(aux.edge_to_qubit) + 1)
        end
    end

    iter = 0
    while true
        basis = cycle_basis(g)
        isempty(basis) && break
        n, idx_longest = findmax(length, basis)
        cycle = basis[idx_longest]
        n > max_cycle_len || break

        iter += 1
        iter ≤ max_iterations ||
            error("cellulate_long_cycles!: max_iterations=$max_iterations " *
                  "reached; longest remaining cycle length = $n")

        u = cycle[1]
        v = cycle[(n ÷ 2) + 1]
        u == v && error(
            "cellulate_long_cycles!: midpoint chord is a self-loop; " *
            "longest cycle = $cycle has malformed length $n")
        if has_edge(g, u, v)
            error("cellulate_long_cycles!: midpoint chord ($u, $v) " *
                  "already exists; longest cycle = $cycle")
        end

        add_edge!(g, u, v)
        push!(aux.edge_to_qubit, length(aux.edge_to_qubit) + 1)
    end

    final_basis = cycle_basis(g)
    m_edges = ne(g)
    edge_index = Dict{Tuple{Int,Int},Int}()
    for (i, e) in enumerate(edges(g))
        edge_index[(min(src(e), dst(e)), max(src(e), dst(e)))] = i
    end
    I = Int[]; J = Int[]
    for (c, cycle) in enumerate(final_basis)
        n_c = length(cycle)
        for i in 1:n_c
            u = cycle[i]; v = cycle[mod1(i + 1, n_c)]
            key = (min(u, v), max(u, v))
            haskey(edge_index, key) ||
                error("cellulate_long_cycles!: edge ($u,$v) of basis cycle $c " *
                      "not found in graph; basis inconsistent with graph")
            push!(I, c); push!(J, edge_index[key])
        end
    end
    n_cycles = length(final_basis)
    aux_cb = sparse(I, J, trues(length(I)), n_cycles, m_edges)
    return AuxiliaryGraph(
        aux.graph,
        aux.port_function,
        aux.edge_to_qubit,
        aux_cb,
        aux.logical_support,
        aux.stabilizer_matchings,
    )
end

"""
Exact relative expansion `β_t(G, U)` (paper Definition 2) by brute-force
enumeration of `2^n` vertex subsets. Returns `typemax(Int) // 1` as +∞ when
the constraint is vacuous. Errors for `n > 24`.
"""
function relative_expansion(graph::SimpleGraph{Int},
                             port_vertices::Vector{Int},
                             t::Int)::Rational{Int}
    n = nv(graph)
    n ≤ 24 || throw(ArgumentError(
        "relative_expansion: brute force only supports n ≤ 24 (got n=$n); " *
        "use sampling-based bounds for larger graphs"))
    t ≥ 1 || throw(ArgumentError("relative_expansion: t must be ≥ 1 (got $t)"))
    U = length(port_vertices)
    port_mask = falses(n)
    for q in port_vertices
        1 ≤ q ≤ n ||
            throw(ArgumentError("port vertex $q out of range 1:$n"))
        port_mask[q] = true
    end
    edge_list = Tuple{Int,Int}[]
    for e in edges(graph)
        push!(edge_list, (src(e), dst(e)))
    end
    best = typemax(Int) // 1
    v_set = falses(n)
    upper = (1 << n)            # 2^n
    @inbounds for k in 0:(upper - 1)
        @inbounds for i in 1:n
            v_set[i] = ((k >> (i - 1)) & 1) == 1
        end
        boundary = 0
        @inbounds for (a, b) in edge_list
            if v_set[a] != v_set[b]
                boundary += 1
            end
        end
        u_size = 0
        @inbounds for i in 1:n
            if port_mask[i] && v_set[i]
                u_size += 1
            end
        end
        denom = min(t, u_size, U - u_size)
        denom > 0 || continue
        ratio = boundary // denom
        if ratio < best
            best = ratio
        end
    end
    best
end

"""
Check `β_d(G, im f) ≥ 1` (paper Theorem 5, desideratum 4). `im f` is the
union of port-function vertices. If the graph fails the bound, the paper
prescribes adding expansion-boosting edges; that fixup is not implemented
since the worked examples satisfy the bound without modification.
"""
function verify_desideratum_4(aux::AuxiliaryGraph, d::Int)::Bool
    pv_set = Set{Int}()
    for vs in aux.port_function
        for v in vs
            push!(pv_set, v)
        end
    end
    port_vertices = sort(collect(pv_set))
    β = relative_expansion(aux.graph, port_vertices, d)
    β ≥ 1
end
