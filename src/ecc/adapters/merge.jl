"""
    incidence_matrix_transpose(aux::AuxiliaryGraph) :: SparseMatrixCSC{Bool, Int}

`(nv, ne)` incidence-matrix transpose of `aux.graph`: row `v`, column `e`
is `true` iff `v` is an endpoint of edge `e`. Edge ordering matches
`edges(aux.graph)`. This is the `G^⊤` block of Figure 7
[swaroop2026universal](@cite), used in the V_l / V_r vertex-Z rows.
"""
function incidence_matrix_transpose(aux::AuxiliaryGraph)::SparseMatrixCSC{Bool, Int}
    n = nv(aux.graph)
    m = ne(aux.graph)
    I = Int[]; J = Int[]
    sizehint!(I, 2 * m); sizehint!(J, 2 * m)
    @inbounds for (e_idx, edge) in enumerate(edges(aux.graph))
        u, v = src(edge), dst(edge)
        push!(I, u); push!(J, e_idx)
        push!(I, v); push!(J, e_idx)
    end
    sparse(I, J, trues(length(I)), n, m)
end

"""
    stabilizer_modification_matrix(aux::AuxiliaryGraph, n_stabilizers::Int) :: SparseMatrixCSC{Bool, Int}

Stabilizer-modification matrix `M` of Figure 7: `M[s, e] = 1` iff edge `e`
is in μ(L_s). Shape `(n_stabilizers, ne(aux.graph))`. Resolves vertex pairs
in `aux.stabilizer_matchings` to current edge indices via a precomputed
lookup, so it works both pre- and post-cellulation. Errors if a stored
pair is not currently an edge (indicates external mutation of the graph).
"""
function stabilizer_modification_matrix(aux::AuxiliaryGraph,
                                          n_stabilizers::Int)::SparseMatrixCSC{Bool, Int}
    m = ne(aux.graph)
    length(aux.stabilizer_matchings) == n_stabilizers ||
        throw(DimensionMismatch(
            "stabilizer_modification_matrix: aux.stabilizer_matchings has " *
            "length $(length(aux.stabilizer_matchings)) but n_stabilizers=$n_stabilizers"))
    # Precomputed pair → edge index — O(m) build vs O(m·matchings) for per-pair lookup.
    pair_idx = Dict{Tuple{Int,Int},Int}()
    @inbounds for (i, e) in enumerate(edges(aux.graph))
        u, v = minmax(src(e), dst(e))
        pair_idx[(u, v)] = i
    end
    I = Int[]; J = Int[]
    @inbounds for s in 1:n_stabilizers
        for pair in aux.stabilizer_matchings[s]
            haskey(pair_idx, pair) ||
                error("stabilizer_modification_matrix: pair $pair for " *
                      "stabilizer $s not found in current graph; " *
                      "AuxiliaryGraph may have been mutated externally")
            push!(I, s); push!(J, pair_idx[pair])
        end
    end
    sparse(I, J, trues(length(I)), n_stabilizers, m)
end

"""
    canonical_H_R(n::Int) :: SparseMatrixCSC{Bool, Int}

The `(n-1, n)` banded full-rank repetition parity-check matrix with
`H_R[i, i] = H_R[i, i+1] = 1`. This is the bridge × adapter block of the
merged H_X (Eq. 37 of [swaroop2026universal](@cite)) and the algorithm-2
target of [`skiptree`](@ref).
"""
function canonical_H_R(n::Int)::SparseMatrixCSC{Bool, Int}
    n ≥ 2 || throw(ArgumentError("canonical_H_R: n must be ≥ 2 (got $n)"))
    rows = n - 1
    I = Vector{Int}(undef, 2 * rows)
    J = Vector{Int}(undef, 2 * rows)
    @inbounds for i in 1:rows
        I[2i - 1] = i; J[2i - 1] = i
        I[2i]     = i; J[2i]     = i + 1
    end
    sparse(I, J, trues(2 * rows), rows, n)
end

"""
    adapter_bridge_x_check_matrix(
        skiptree_l::SkipTreeOutput,
        skiptree_r::SkipTreeOutput,
        adapter_width::Int,
    ) :: SparseMatrixCSC{Bool, Int}

Bridge X-check matrix ``[T_l \\mid H_R(w) \\mid T_r]`` (Eq. 37 of
[swaroop2026universal](@cite)). Shape `(w - 1, m_l + w + m_r)`. Errors with
`DimensionMismatch` if `T_l` or `T_r` does not have `w - 1` rows.
"""
function adapter_bridge_x_check_matrix(
    skiptree_l::SkipTreeOutput,
    skiptree_r::SkipTreeOutput,
    adapter_width::Int,
)::SparseMatrixCSC{Bool, Int}
    adapter_width ≥ 2 ||
        throw(ArgumentError("adapter_bridge_x_check_matrix: adapter_width must be ≥ 2 (got $adapter_width)"))
    expected_rows = adapter_width - 1
    size(skiptree_l.T, 1) == expected_rows || throw(DimensionMismatch(
        "adapter_bridge_x_check_matrix: skiptree_l.T has $(size(skiptree_l.T,1)) " *
        "rows; expected $expected_rows (= adapter_width - 1)"))
    size(skiptree_r.T, 1) == expected_rows || throw(DimensionMismatch(
        "adapter_bridge_x_check_matrix: skiptree_r.T has $(size(skiptree_r.T,1)) " *
        "rows; expected $expected_rows (= adapter_width - 1)"))
    hcat(skiptree_l.T, canonical_H_R(adapter_width), skiptree_r.T)
end

"""
    assemble_merged_code_intercode(
        code1::CSS, code2::CSS,
        aux1::AuxiliaryGraph, aux2::AuxiliaryGraph,
        skiptree1::SkipTreeOutput, skiptree2::SkipTreeOutput,
    ) :: CSS

Inter-code merged CSS code (Figure 7 / Eq. 37 of [swaroop2026universal](@cite));
`code1` and `code2` must be distinct objects. The two `aux.logical_support`
must have equal length (Lemma 9). Auxiliary graphs must have been cellulated
(so `cycle_basis_matrix` is populated, though zero rows is fine for trees).

# Column layout (n = n1 + m1 + w + m2 + n2)

| code1 | G_1 edges | A adapter | G_2 edges | code2 |
|-------|-----------|-----------|-----------|-------|
| 1..n1 | n1+1..n1+m1 | next w | next m2 | last n2 |

# H_X row blocks

| Row block       | code1 | G_1 | A      | G_2 | code2 |
|-----------------|-------|-----|--------|-----|-------|
| code1 X-stab    | Hx_1  | M_l | 0      | 0   | 0     |
| code2 X-stab    | 0     | 0   | 0      | M_r | Hx_2  |
| N_l cycles      | 0     | N_l | 0      | 0   | 0     |
| bridge (w-1)    | 0     | T_l | H_R(w) | T_r | 0     |
| N_r cycles      | 0     | 0   | 0      | N_r | 0     |

# H_Z row blocks

| Row block      | code1 | G_1   | A   | G_2   | code2 |
|----------------|-------|-------|-----|-------|-------|
| code1 Z-stab   | Hz_1  | 0     | 0   | 0     | 0     |
| code2 Z-stab   | 0     | 0     | 0   | 0     | Hz_2  |
| V_l vertex (w) | F_l   | G_l^T | P_l | 0     | 0     |
| V_r vertex (w) | 0     | 0     | P_r | G_r^T | F_r   |

The `P_l`/`P_r` A-blocks are the SkipTree permutations and are required
for CSS commutation between bridge X-checks and V-vertex Z-checks.
"""
function assemble_merged_code_intercode(
    code1::CSS, code2::CSS,
    aux1::AuxiliaryGraph, aux2::AuxiliaryGraph,
    skiptree1::SkipTreeOutput, skiptree2::SkipTreeOutput,
)::CSS
    w_l = length(aux1.logical_support)
    w_r = length(aux2.logical_support)
    w_l == w_r ||
        throw(DimensionMismatch(
            "assemble_merged_code_intercode: aux1.logical_support has " *
            "length $w_l, aux2.logical_support has length $w_r; " *
            "adapter requires equal port widths"))
    w = w_l
    # Empty cycle_basis_matrix (tree-like aux) is valid: just gives zero N_l/N_r rows.
    nv(aux1.graph) == w || throw(DimensionMismatch(
        "aux1 has $(nv(aux1.graph)) vertices; expected $w = |logical_support|"))
    nv(aux2.graph) == w || throw(DimensionMismatch(
        "aux2 has $(nv(aux2.graph)) vertices; expected $w = |logical_support|"))

    n1 = size(code1.Hx, 2)
    n2 = size(code2.Hx, 2)
    m1 = ne(aux1.graph)
    m2 = ne(aux2.graph)
    sx1 = size(code1.Hx, 1); sz1 = size(code1.Hz, 1)
    sx2 = size(code2.Hx, 1); sz2 = size(code2.Hz, 1)

    M_l  = stabilizer_modification_matrix(aux1, sx1)
    M_r  = stabilizer_modification_matrix(aux2, sx2)
    N_l  = aux1.cycle_basis_matrix
    N_r  = aux2.cycle_basis_matrix
    HR   = canonical_H_R(w)
    T_l  = skiptree1.T
    T_r  = skiptree2.T
    F_l  = port_function_as_matrix(aux1, n1)
    F_r  = port_function_as_matrix(aux2, n2)
    GT_l = incidence_matrix_transpose(aux1)
    GT_r = incidence_matrix_transpose(aux2)
    P_l  = skiptree1.P
    P_r  = skiptree2.P

    Hx_BB = sparse(code1.Hx)
    Hz_BB = sparse(code1.Hz)
    Hx_LP = sparse(code2.Hx)
    Hz_LP = sparse(code2.Hz)

    Z = (r::Int, c::Int) -> spzeros(Bool, r, c)
    nc_l = size(N_l, 1)
    nc_r = size(N_r, 1)

    Hx_row1 = hcat(Hx_BB,        M_l,          Z(sx1, w),    Z(sx1, m2),  Z(sx1, n2))
    Hx_row2 = hcat(Z(sx2, n1),   Z(sx2, m1),   Z(sx2, w),    M_r,         Hx_LP)
    Hx_row3 = hcat(Z(nc_l, n1),  N_l,          Z(nc_l, w),   Z(nc_l, m2), Z(nc_l, n2))
    Hx_row4 = hcat(Z(w - 1, n1), T_l,          HR,           T_r,         Z(w - 1, n2))
    Hx_row5 = hcat(Z(nc_r, n1),  Z(nc_r, m1),  Z(nc_r, w),   N_r,         Z(nc_r, n2))
    Hx_merged = vcat(Hx_row1, Hx_row2, Hx_row3, Hx_row4, Hx_row5)

    Hz_row1 = hcat(Hz_BB,        Z(sz1, m1),  Z(sz1, w),   Z(sz1, m2),  Z(sz1, n2))
    Hz_row2 = hcat(Z(sz2, n1),   Z(sz2, m1),  Z(sz2, w),   Z(sz2, m2),  Hz_LP)
    Hz_row3 = hcat(F_l,          GT_l,        P_l,         Z(w, m2),    Z(w, n2))
    Hz_row4 = hcat(Z(w, n1),     Z(w, m1),    P_r,         GT_r,        F_r)
    Hz_merged = vcat(Hz_row1, Hz_row2, Hz_row3, Hz_row4)

    CSS(Matrix(Hx_merged), Matrix(Hz_merged))
end

"""
    build_adapter_intercode(
        pair::CodePair;
        max_cycle_len_1::Union{Int,Nothing} = nothing,
        max_cycle_len_2::Union{Int,Nothing} = nothing,
        chords_1::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
        chords_2::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
    ) :: Adapter

Inter-code joint Z̄_1 Z̄_2 measurement. Builds both aux graphs,
cellulates, runs SkipTree, and assembles per
[`assemble_merged_code_intercode`](@ref).

When `max_cycle_len_i = nothing`, defaults to the max X-stabilizer weight
of the corresponding code — matches paper §VII.A's heuristic ("LP_2 has
weight-7 stabilizers, so we choose not to cellulate") and reproduces
Example B (BB[[98,6,12]] × LP[[200,20,10]]) with no kwargs (`n=355`).
`chords_{1,2}` forward to [`cellulate_long_cycles!`](@ref) when you need
byte-identical reproduction.

Errors `ArgumentError` if `code1 === code2` (use the intra-code path),
`DimensionMismatch` if the two logical supports differ in length.
"""
function build_adapter_intercode(
    pair::CodePair;
    max_cycle_len_1::Union{Int,Nothing} = nothing,
    max_cycle_len_2::Union{Int,Nothing} = nothing,
    chords_1::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
    chords_2::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
)::Adapter
    pair.code1 !== pair.code2 ||
        throw(ArgumentError("build_adapter_intercode: code1 and code2 must be distinct objects; use intra-code path for measurements on the same codeblock"))
    w = length(pair.logical1_support)
    length(pair.logical2_support) == w || throw(DimensionMismatch(
        "build_adapter_intercode: logical supports have unequal lengths " *
        "($w vs $(length(pair.logical2_support)))"))
    mcl1 = max_cycle_len_1 === nothing ?
           Int(maximum(sum(pair.code1.Hx; dims = 2))) : max_cycle_len_1
    mcl2 = max_cycle_len_2 === nothing ?
           Int(maximum(sum(pair.code2.Hx; dims = 2))) : max_cycle_len_2
    aux1 = cellulate_long_cycles!(
        build_initial_aux_graph(pair.logical1_support, pair.code1);
        max_cycle_len = mcl1, chords = chords_1)
    aux2 = cellulate_long_cycles!(
        build_initial_aux_graph(pair.logical2_support, pair.code2);
        max_cycle_len = mcl2, chords = chords_2)
    st1 = skiptree(aux1.graph)
    st2 = skiptree(aux2.graph)
    merged = assemble_merged_code_intercode(pair.code1, pair.code2,
                                             aux1, aux2, st1, st2)
    Adapter(pair, aux1, aux2, st1, st2, w, merged)
end

"""
    assemble_merged_code_intracode(
        code::CSS,
        aux1::AuxiliaryGraph, aux2::AuxiliaryGraph,
        skiptree1::SkipTreeOutput, skiptree2::SkipTreeOutput,
    ) :: CSS

Intra-code merged CSS code (both logicals on the same block) per
Figure 7 / Eq. 37 of [swaroop2026universal](@cite).

# Column layout (n = n_code + m_l + w + m_r)

| code | G_1 edges | A adapter | G_2 edges |
|------|-----------|-----------|-----------|
| 1..n_code | next m_l | next w | last m_r |

where `w = min(|aux1.logical_support|, |aux2.logical_support|)`.

# H_X row blocks

| Row block      | code | G_1 | A      | G_2 |
|----------------|------|-----|--------|-----|
| code X-stab    | H_X  | M_l | 0      | M_r |
| N_l cycles     | 0    | N_l | 0      | 0   |
| bridge (w-1)   | 0    | T_l | H_R(w) | T_r |
| N_r cycles     | 0    | 0   | 0      | N_r |

# H_Z row blocks

| Row block        | code | G_1   | A    | G_2   |
|------------------|------|-------|------|-------|
| code Z-stab      | H_Z  | 0     | 0    | 0     |
| V_l vertex (w_l) | F_l  | G_l^T | P_l' | 0     |
| V_r vertex (w_r) | F_r  | 0     | P_r' | G_r^T |

# Differences vs. inter-code

1. A single code X-stab can touch both logicals, so M_l and M_r share
   the same code-X-stab row block (not two separate blocks).
2. When `w_l ≠ w_r`, the adapter width is `w = min(w_l, w_r)`. We trim
   `T_{l,r}` to `w − 1` rows and `P_{l,r}` to `w` columns; the `w_l − w`
   "extra" aux1 vertices keep their F/G^T entries (needed for vertex-Z
   commutation) but get a zero A-block, and the bridge only touches the
   first `w` label positions so commutation still holds.
3. Both F_l and F_r target the code column block — no separate code2 block.
"""
function assemble_merged_code_intracode(
    code::CSS,
    aux1::AuxiliaryGraph, aux2::AuxiliaryGraph,
    skiptree1::SkipTreeOutput, skiptree2::SkipTreeOutput,
)::CSS
    w_l = length(aux1.logical_support)
    w_r = length(aux2.logical_support)
    nv(aux1.graph) == w_l || throw(DimensionMismatch(
        "aux1 has $(nv(aux1.graph)) vertices; expected $w_l = |logical_support|"))
    nv(aux2.graph) == w_r || throw(DimensionMismatch(
        "aux2 has $(nv(aux2.graph)) vertices; expected $w_r = |logical_support|"))
    w = min(w_l, w_r)
    w ≥ 2 || throw(ArgumentError(
        "assemble_merged_code_intracode: adapter width = $w < 2"))

    n_code = size(code.Hx, 2)
    sx = size(code.Hx, 1)
    sz = size(code.Hz, 1)
    m1 = ne(aux1.graph)
    m2 = ne(aux2.graph)

    M_l = stabilizer_modification_matrix(aux1, sx)
    M_r = stabilizer_modification_matrix(aux2, sx)
    N_l = aux1.cycle_basis_matrix
    N_r = aux2.cycle_basis_matrix
    HR  = canonical_H_R(w)
    # T trimmed to (w-1) rows, P to w columns when w < w_{l,r} (see header).
    T_l = skiptree1.T[1:(w - 1), :]
    T_r = skiptree2.T[1:(w - 1), :]
    F_l  = port_function_as_matrix(aux1, n_code)
    F_r  = port_function_as_matrix(aux2, n_code)
    GT_l = incidence_matrix_transpose(aux1)
    GT_r = incidence_matrix_transpose(aux2)
    P_l = skiptree1.P[:, 1:w]
    P_r = skiptree2.P[:, 1:w]

    Hx_c = sparse(code.Hx)
    Hz_c = sparse(code.Hz)

    Z = (r::Int, c::Int) -> spzeros(Bool, r, c)
    nc_l = size(N_l, 1)
    nc_r = size(N_r, 1)

    Hx_row1 = hcat(Hx_c,            M_l,         Z(sx, w),    M_r)             # code X-stab
    Hx_row2 = hcat(Z(nc_l, n_code), N_l,         Z(nc_l, w),  Z(nc_l, m2))     # N_l
    Hx_row3 = hcat(Z(w-1,  n_code), T_l,         HR,          T_r)             # bridge
    Hx_row4 = hcat(Z(nc_r, n_code), Z(nc_r, m1), Z(nc_r, w),  N_r)             # N_r
    Hx_merged = vcat(Hx_row1, Hx_row2, Hx_row3, Hx_row4)

    Hz_row1 = hcat(Hz_c,    Z(sz, m1),  Z(sz, w),  Z(sz, m2))        # code Z-stab
    Hz_row2 = hcat(F_l,     GT_l,       P_l,       Z(w_l, m2))       # V_l (w_l rows)
    Hz_row3 = hcat(F_r,     Z(w_r, m1), P_r,       GT_r)             # V_r (w_r rows)
    Hz_merged = vcat(Hz_row1, Hz_row2, Hz_row3)

    CSS(Matrix(Hx_merged), Matrix(Hz_merged))
end

"""
    build_adapter_intracode(
        pair::CodePair;
        max_cycle_len_1::Union{Int,Nothing} = nothing,
        max_cycle_len_2::Union{Int,Nothing} = nothing,
        chords_1::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
        chords_2::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
    ) :: Adapter

Intra-code (same codeblock) version of [`build_adapter_intercode`](@ref).
Reproduces paper Example C ([[150, 5, 12]], BB intra with Z_1/Z_3, adapter
width 12) with default kwargs. Errors `ArgumentError` if
`code1 !== code2` or if either support has < 2 qubits.
"""
function build_adapter_intracode(
    pair::CodePair;
    max_cycle_len_1::Union{Int,Nothing} = nothing,
    max_cycle_len_2::Union{Int,Nothing} = nothing,
    chords_1::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
    chords_2::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
)::Adapter
    pair.code1 === pair.code2 ||
        throw(ArgumentError("build_adapter_intracode: code1 and code2 must be the same object; use build_adapter_intercode for distinct codeblocks"))
    code = pair.code1
    w = min(length(pair.logical1_support), length(pair.logical2_support))
    w ≥ 2 || throw(ArgumentError(
        "build_adapter_intracode: adapter width = $w < 2"))
    mcl1 = max_cycle_len_1 === nothing ?
           Int(maximum(sum(code.Hx; dims = 2))) : max_cycle_len_1
    mcl2 = max_cycle_len_2 === nothing ?
           Int(maximum(sum(code.Hx; dims = 2))) : max_cycle_len_2
    aux1 = cellulate_long_cycles!(
        build_initial_aux_graph(pair.logical1_support, code);
        max_cycle_len = mcl1, chords = chords_1)
    aux2 = cellulate_long_cycles!(
        build_initial_aux_graph(pair.logical2_support, code);
        max_cycle_len = mcl2, chords = chords_2)
    st1 = skiptree(aux1.graph)
    st2 = skiptree(aux2.graph)
    merged = assemble_merged_code_intracode(code, aux1, aux2, st1, st2)
    Adapter(pair, aux1, aux2, st1, st2, w, merged)
end

"""
    build_adapter(pair::CodePair; kwargs...) :: Adapter

Dispatches to [`build_adapter_intracode`](@ref) if `code1 === code2`,
else [`build_adapter_intercode`](@ref). Reproduces paper Examples B and
C with no extra kwargs. Use [`AdapterMergedCode`](@ref) to wrap the
result for the standard `AbstractECC` dispatch.

```jldoctest adapter_examples
julia> using QuantumClifford, QECCore

julia> using QuantumClifford.ECC: Surface, logz_ops, code_n, code_k

julia> using QuantumClifford.ECC: parity_matrix_x, parity_matrix_z

julia> using QuantumClifford: stab_to_gf2

julia> c1 = CSS(Matrix{Bool}(parity_matrix_x(Surface(3,3))), Matrix{Bool}(parity_matrix_z(Surface(3,3))));

julia> c2 = CSS(Matrix{Bool}(parity_matrix_x(Surface(3,3))), Matrix{Bool}(parity_matrix_z(Surface(3,3))));

julia> z = sort(findall(!iszero, stab_to_gf2(logz_ops(Surface(3,3)))[1, 14:26]));

julia> adapter = build_adapter(CodePair(c1, c2, z, z));

julia> (code_n(adapter.merged_code), code_k(adapter.merged_code), adapter.adapter_width)
(33, 1, 3)
```
"""
function build_adapter(pair::CodePair; kwargs...)::Adapter
    if pair.code1 === pair.code2
        build_adapter_intracode(pair; kwargs...)
    else
        build_adapter_intercode(pair; kwargs...)
    end
end

# AbstractCSSCode interface for AdapterMergedCode — delegates to the inner CSS.
parity_matrix_x(c::AdapterMergedCode) = c.adapter.merged_code.Hx
parity_matrix_z(c::AdapterMergedCode) = c.adapter.merged_code.Hz
code_n(c::AdapterMergedCode) = code_n(c.adapter.merged_code)
code_s(c::AdapterMergedCode) = code_s(c.adapter.merged_code)

"""
    deform_code(
        code::CSS,
        logical_support::Vector{Int};
        max_cycle_len::Union{Int,Nothing} = nothing,
        chords::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
    ) :: CSS

Single-logical deformed code (Figure 3 / Eq. 9 of [swaroop2026universal](@cite)):
attach one auxiliary graph + cycle checks to `code` so the Z operator on
`logical_support` becomes a stabilizer product. No bridge, no SkipTree —
use [`build_adapter`](@ref) for the joint-measurement case.

```
H_X_deform = [ H_X | M  ]      original X-stabs deformed by matchings
             [ 0   | N  ]      cycle X-checks (one per basis cycle)

H_Z_deform = [ H_Z | 0   ]     original Z-stabs unchanged
             [ F   | G^T ]     vertex Z-checks (one per port vertex)
```

`max_cycle_len = nothing` defaults to the original code's max X-stab weight
(matching paper §VII.A's heuristic). `chords` forwards to
[`cellulate_long_cycles!`](@ref) for byte-identical paper reproduction.

# Reproducibility — paper Table II

| Original code | Logical | Deformed params | Default cellulation? |
|---|---|---|:---:|
| BB[[98,6,12]] | Z_1 | [[121,5,12]] | ✓ |
| LP[[200,20,10]] | Z_2 | [[220,19,10]] | ✓ |
| BB[[98,6,12]] (full-rank basis) | Z_3 | [[115,5,12]] | needs `chords=[(1,11)]` |

For BB / Z_3, our default midpoint picks chord `(4, 9)`; Zenodo's
NetworkX pipeline picks `(1, 11)`. Both give valid 17-edge cellulations,
but our default produces `d ≤ 11` and Zenodo's gives `d = 12`.
Single-logical reproduction therefore needs the explicit `chords` kwarg.
The joint-measurement case via [`build_adapter`](@ref) is unaffected —
the bridge X-checks absorb the extra weight-`(d-1)` logical that the
chord choice introduces.

```jldoctest
julia> using QuantumClifford, QECCore

julia> using QuantumClifford.ECC: Surface, logz_ops, code_n, code_k

julia> using QuantumClifford.ECC: parity_matrix_x, parity_matrix_z

julia> using QuantumClifford: stab_to_gf2

julia> c = CSS(Matrix{Bool}(parity_matrix_x(Surface(3,3))), Matrix{Bool}(parity_matrix_z(Surface(3,3))));

julia> z = sort(findall(!iszero, stab_to_gf2(logz_ops(Surface(3,3)))[1, 14:26]));

julia> dc = deform_code(c, z);

julia> (code_n(dc), code_k(dc))
(15, 0)
```
"""
function deform_code(
    code::CSS,
    logical_support::Vector{Int};
    max_cycle_len::Union{Int,Nothing} = nothing,
    chords::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
)::CSS
    mcl = max_cycle_len === nothing ?
          Int(maximum(sum(code.Hx; dims = 2))) : max_cycle_len
    aux = cellulate_long_cycles!(
        build_initial_aux_graph(logical_support, code);
        max_cycle_len = mcl, chords = chords)

    n_code = size(code.Hx, 2)
    sx = size(code.Hx, 1)
    sz = size(code.Hz, 1)
    m  = ne(aux.graph)

    M  = stabilizer_modification_matrix(aux, sx)
    N  = aux.cycle_basis_matrix
    F  = port_function_as_matrix(aux, n_code)
    GT = incidence_matrix_transpose(aux)
    n_c = size(N, 1)

    Hx_code = sparse(code.Hx)
    Hz_code = sparse(code.Hz)
    Z = (r::Int, c::Int) -> spzeros(Bool, r, c)

    Hx_row1 = hcat(Hx_code,        M)
    Hx_row2 = hcat(Z(n_c, n_code), N)
    Hx_deform = vcat(Hx_row1, Hx_row2)

    Hz_row1 = hcat(Hz_code, Z(sz, m))
    Hz_row2 = hcat(F,       GT)
    Hz_deform = vcat(Hz_row1, Hz_row2)

    CSS(Matrix(Hx_deform), Matrix(Hz_deform))
end
