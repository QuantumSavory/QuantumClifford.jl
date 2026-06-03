"""
    build_initial_aux_graph(logical_support::Vector{Int}, code::CSS) :: AuxiliaryGraph

Construct the initial auxiliary graph G_0 (§II.C of [swaroop2026universal](@cite))
for measuring Z̄ with the given `logical_support` on `code`. One vertex per
support qubit (identity port function), one edge per pair of consecutive
vertices in `sort(supp(s) ∩ logical_support)` for each X-stabilizer row `s`
that overlaps the logical. If the resulting graph is disconnected, emits a
`@warn` and greedily adds repair edges between components.

Satisfies desiderata 0–2 of Theorem 5 (connectivity, bounded degree, short
matchings). Desideratum 3 (sparse cycle basis) needs
[`cellulate_long_cycles!`](@ref); desideratum 4 (relative expansion ≥ 1)
must be checked separately and only requires a fixup if it fails — the
worked examples in paper §VII already meet it without modification.

`edge_to_qubit` is filled with `1:ne(graph)` as a placeholder; the real
global-qubit assignment happens at merge time. `cycle_basis_matrix` is
empty until cellulation runs.

The sorted-adjacent-pair matching is deterministic and reproduces Zenodo
on the small worked examples, but is not optimal in general; a smarter
matching that minimises edge reuse across stabilizers may give better
expansion downstream.
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
    port_function_as_matrix(aux::AuxiliaryGraph, n_qubits::Int) :: SparseMatrixCSC{Bool, Int}

Return matrix `F` of Figure 7 of [swaroop2026universal](@cite):
`F[v, q] = 1` iff `v ∈ port_function[q_idx]` where `q_idx` is `q`'s position
in `aux.logical_support`. Shape `(nv(aux.graph), n_qubits)`.
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
    edge_pair_to_index(graph::SimpleGraph{Int}, pair::Tuple{Int,Int}) :: Int

Return the position of the edge `pair = (u, v)` (with `u ≤ v`) in
`edges(graph)` iteration order; returns `0` if the pair is not an edge.
This iteration order is the canonical column convention for the merged
H_X / H_Z edge blocks.

Linear scan, O(`ne(graph)`). Callers in hot loops should precompute a
`Dict{Tuple{Int,Int}, Int}` once and reuse it (see
[`stabilizer_modification_matrix`](@ref)).
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
    cellulate_long_cycles!(aux::AuxiliaryGraph; max_cycle_len=6,
                            max_iterations=100, chords=nothing) :: AuxiliaryGraph

Cellulate `aux.graph` so every basis cycle has length ≤ `max_cycle_len`
(desideratum 3a of [swaroop2026universal](@cite)). Default algorithm:
repeatedly take the longest cycle in `Graphs.cycle_basis` and add the
chord `(cycle[1], cycle[n÷2 + 1])` (paper §VII.A midpoint split). If
`chords` is supplied, those are applied in order first, then the
midpoint algorithm runs on whatever remains.

Mutates `aux.graph`, extends `aux.edge_to_qubit` with placeholders for
each new edge, and rebuilds `aux.cycle_basis_matrix` from the final
basis. Returns a new `AuxiliaryGraph` sharing the mutated graph (the
struct is immutable except for `graph` and `edge_to_qubit`, which is
why we can't rebind `cycle_basis_matrix` in place).

# Errors

`ArgumentError` for out-of-range chord endpoints, self-loops, or chord
already present. `error` if `max_iterations` is exhausted (likely a
degenerate input).

# Why the `chords` kwarg exists — distance sensitivity to basis choice

`Graphs.cycle_basis` and NetworkX choose different DFS bases for the
same graph, so the "longest cycle" we split differs from Zenodo's.
Only the final edge count, the basis length bound, `N · G ≡ 0 (mod 2)`,
and the rank `m − n + 1` are guaranteed.

This matters: paper Table II's `BB / Z_3` deformed code is `[[115, 5, 12]]`
with Zenodo's chord `(1, 11)` but `[[115, 5, ≤11]]` with our default
midpoint chord `(4, 9)`. For byte-identical paper reproduction of the
single-logical case, pass `chords=[(1, 11)]`. For joint measurements via
[`build_adapter`](@ref), the bridge X-checks absorb the extra weight-
`(d-1)` logical empirically, so the merged code's distance comes out
right regardless.
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
    relative_expansion(graph::SimpleGraph{Int}, port_vertices::Vector{Int}, t::Int) :: Rational{Int}

Exact relative expansion `β_t(G, U)` (Definition 2 of
[swaroop2026universal](@cite)) by brute-force over `2^n` vertex subsets.
Returns `typemax(Int) // 1` as +∞ when no subset gives a positive
denominator. Errors for `n > 24` — use an LP relaxation or sampling
bound for larger graphs.
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
    verify_desideratum_4(aux::AuxiliaryGraph, d::Int) :: Bool

Check `β_d(G, im f) ≥ 1` (desideratum 4 of Theorem 5,
[swaroop2026universal](@cite)). `im f` is the union of port-function
vertices. The paper's §VII.A worked examples all satisfy this without
modification. Graphs that don't need expansion-boosting edges per §II.C
(not implemented yet).
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
