"""
The `G^⊤` block from Figure 7 [swaroop2026universal](@cite). Row `v`, column
`e` is true iff `v` is an endpoint of edge `e`, in `edges(aux.graph)` order.
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
The `M` block from Figure 7 [swaroop2026universal](@cite). `M[s, e] = 1` iff
edge `e` is in the matching μ(L_s) of X-stabilizer `s`. Stored vertex pairs
are resolved to current edge indices at call time, so it stays valid across
cellulation chord additions.
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
The `(n-1, n)` banded full-rank repetition parity-check matrix:
`H_R[i, i] = H_R[i, i+1] = 1`. Bridge × adapter block of the merged H_X
(paper Eq. 37) and the target of [`skiptree`](@ref).
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
The bridge X-check block ``[T_l \\mid H_R(w) \\mid T_r]`` from paper Eq. 37.
Shape `(w - 1, m_l + w + m_r)`.
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
Inter-code merged CSS code (paper Figure 7 / Eq. 37). `code1` and `code2` must
be distinct objects with equal-length logical supports (Lemma 9). Auxiliary
graphs must have been cellulated.

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

The `P_l`/`P_r` A-blocks are the SkipTree permutations, required for CSS
commutation between bridge X-checks and V-vertex Z-checks.
"""
function assemble_merged_code_intercode(
    code1, code2,
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

    Hx_1 = Matrix{Bool}(parity_matrix_x(code1))
    Hz_1 = Matrix{Bool}(parity_matrix_z(code1))
    Hx_2 = Matrix{Bool}(parity_matrix_x(code2))
    Hz_2 = Matrix{Bool}(parity_matrix_z(code2))

    n1 = size(Hx_1, 2)
    n2 = size(Hx_2, 2)
    m1 = ne(aux1.graph)
    m2 = ne(aux2.graph)
    sx1 = size(Hx_1, 1); sz1 = size(Hz_1, 1)
    sx2 = size(Hx_2, 1); sz2 = size(Hz_2, 1)

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

    Hx_BB = sparse(Hx_1)
    Hz_BB = sparse(Hz_1)
    Hx_LP = sparse(Hx_2)
    Hz_LP = sparse(Hz_2)

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
Inter-code joint Z̄_1 Z̄_2 measurement. Builds both aux graphs, cellulates,
runs SkipTree, and assembles via [`assemble_merged_code_intercode`](@ref).

`max_cycle_len_i = nothing` defaults to the max X-stabilizer weight of
the corresponding code (paper §VII.A heuristic). With no kwargs this
reproduces Example B (BB[[98,6,12]] × LP[[200,20,10]] → n=355).
`chords_{1,2}` forwards to [`cellulate_long_cycles!`](@ref).
"""
function _build_adapter_intercode(
    pair::CodePair;
    max_cycle_len_1::Union{Int,Nothing} = nothing,
    max_cycle_len_2::Union{Int,Nothing} = nothing,
    chords_1::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
    chords_2::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
)::Adapter
    w = length(pair.logical1_support)
    length(pair.logical2_support) == w || throw(DimensionMismatch(
        "build_adapter: logical supports have unequal lengths " *
        "($w vs $(length(pair.logical2_support)))"))
    Hx_1 = parity_matrix_x(pair.code1)
    Hx_2 = parity_matrix_x(pair.code2)
    mcl1 = max_cycle_len_1 === nothing ?
           Int(maximum(sum(Hx_1; dims = 2))) : max_cycle_len_1
    mcl2 = max_cycle_len_2 === nothing ?
           Int(maximum(sum(Hx_2; dims = 2))) : max_cycle_len_2
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
Intra-code merged CSS code (both logicals on the same block; paper Figure 7 /
Eq. 37).

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

1. A single code X-stab can touch both logicals, so M_l and M_r share the
   code-X-stab row block.
2. When `w_l ≠ w_r`, the adapter width is `w = min(w_l, w_r)`. `T_{l,r}` is
   trimmed to `w − 1` rows and `P_{l,r}` to `w` columns. The `w_l − w`
   "extra" aux1 vertices keep their F/G^T entries (needed for vertex-Z
   commutation) but get a zero A-block; the bridge only touches the first
   `w` label positions so commutation still holds.
3. Both F_l and F_r target the code column block (no separate code2).
"""
function assemble_merged_code_intracode(
    code,
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

    Hx_code = Matrix{Bool}(parity_matrix_x(code))
    Hz_code = Matrix{Bool}(parity_matrix_z(code))
    n_code = size(Hx_code, 2)
    sx = size(Hx_code, 1)
    sz = size(Hz_code, 1)
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

    Hx_c = sparse(Hx_code)
    Hz_c = sparse(Hz_code)

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
Intra-code joint Z̄_1 Z̄_2 measurement on a single codeblock. Reproduces
paper Example C ([[150, 5, 12]], BB intra with Z_1/Z_3, adapter width 12)
with default kwargs.
"""
function _build_adapter_intracode(
    pair::CodePair;
    max_cycle_len_1::Union{Int,Nothing} = nothing,
    max_cycle_len_2::Union{Int,Nothing} = nothing,
    chords_1::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
    chords_2::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
)::Adapter
    code = pair.code1
    w = min(length(pair.logical1_support), length(pair.logical2_support))
    w ≥ 2 || throw(ArgumentError(
        "build_adapter: adapter width = $w < 2"))
    Hx_code = parity_matrix_x(code)
    mcl1 = max_cycle_len_1 === nothing ?
           Int(maximum(sum(Hx_code; dims = 2))) : max_cycle_len_1
    mcl2 = max_cycle_len_2 === nothing ?
           Int(maximum(sum(Hx_code; dims = 2))) : max_cycle_len_2
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
Build a universal qLDPC adapter between two CSS codeblocks for joint
``\\bar Z_1 \\bar Z_2`` measurement, following [swaroop2026universal](@cite).
The returned [`Adapter`](@ref) is an `AbstractCSSCode`, usable directly
with `code_n`, `parity_matrix_x`, `distance`, and the decoder harnesses.

Two call signatures:

- `build_adapter(code1, code2, logical1, logical2; kwargs...)` — inter-code,
  fuses one logical from each codeblock. `code1 !== code2` is required.
- `build_adapter(code, logical1, logical2; kwargs...)` — intra-code, fuses
  two logicals on the same codeblock.

`logical1` and `logical2` may be passed either as a `Vector{Int}` of qubit
indices, or as a `Z`-only `PauliOperator` (the qubit support is extracted
internally; non-`Z` Paulis are rejected).

`code`, `code1`, `code2` only need to implement `parity_matrix_x` and
`parity_matrix_z`.

Optional kwargs:

- `max_cycle_len_1`, `max_cycle_len_2`: per-side cellulation cycle-length
  budget (paper §VII.A heuristic, default = max X-stabilizer weight).
- `chords_1`, `chords_2`: per-side cellulation chord choices (paper
  Example C `Z_3` reproducibility — see [`cellulate_long_cycles!`](@ref)).

```jldoctest adapter_examples
julia> using QuantumClifford, QECCore

julia> using QuantumClifford.ECC

julia> c1 = Surface(3, 3);

julia> c2 = Surface(3, 3);

julia> z = logz_ops(Surface(3, 3))[1];

julia> adapter = build_adapter(c1, c2, z, z);

julia> (code_n(adapter), code_k(adapter), adapter.adapter_width)
(33, 1, 3)
```
"""
function build_adapter end

# Inter-code: two independent codeblocks.
function build_adapter(code1, code2, logical1, logical2; kwargs...)::Adapter
    s1 = _as_support(logical1, code1, "logical1")
    s2 = _as_support(logical2, code2, "logical2")
    _build_adapter_intercode(CodePair(code1, code2, s1, s2; intracode=false); kwargs...)
end

# Intra-code: single codeblock with two logicals.
function build_adapter(code, logical1, logical2; kwargs...)::Adapter
    s1 = _as_support(logical1, code, "logical1")
    s2 = _as_support(logical2, code, "logical2")
    _build_adapter_intracode(CodePair(code, code, s1, s2; intracode=true); kwargs...)
end

# Convert a logical-operator argument to a sorted qubit-support vector.
_as_support(supp::AbstractVector{<:Integer}, code, label) = collect(Int, supp)

function _as_support(p::PauliOperator, code, label)
    n = size(parity_matrix_x(code), 2)
    nqubits(p) == n || throw(DimensionMismatch(
        "build_adapter: $label has $(nqubits(p)) qubits, code has $n"))
    any(xbit(p)) && throw(ArgumentError(
        "build_adapter: $label must be a Z-only Pauli operator (got non-trivial X support)"))
    sort!(supportz(p))
end

# `Adapter <: AbstractCSSCode` — delegate the CSS interface to the merged code.
parity_matrix_x(a::Adapter) = a.merged_code.Hx
parity_matrix_z(a::Adapter) = a.merged_code.Hz
code_n(a::Adapter) = code_n(a.merged_code)
code_s(a::Adapter) = code_s(a.merged_code)

"""
Row indices into `parity_checks(adapter)` whose product is the joint logical
Pauli being adapted. After standard syndrome extraction on the merged code,
XOR the measurement outcomes at these indices to recover the joint
measurement result.

Inter-code adapter (built via `build_adapter(code1, code2, logical1, logical2)`):
the product equals `Z̄_l ⊗ Z̄_r` embedded in the merged qubit space (`Z̄_l` on
the code1 columns, `Z̄_r` on the code2 columns). Returns the `V_l` and `V_r`
vertex-Z row blocks.

Intra-code adapter (built via `build_adapter(code, logical1, logical2)`): the
product equals `Z̄_l · Z̄_r` (the product of the two Z-logicals on the shared
code block).

The construction: the `V_l` rows sum to `Z̄_l` on the code data qubits
tensored with `Z` on every A-block adapter qubit (the SkipTree permutation
`P_l` puts exactly one `1` per A column); similarly the `V_r` rows sum to
`Z̄_r` tensored with the same A-block `Z`. The A-block contributions cancel
(`Z² = I` over GF(2)), leaving the joint logical on the data qubits.
"""
function joint_logical_recipe(a::Adapter)::Vector{Int}
    n_xstabs = size(parity_matrix_x(a), 1)
    sz_1 = size(parity_matrix_z(a.code_pair.code1), 1)
    if a.code_pair.intracode
        w_l = length(a.code_pair.logical1_support)
        w_r = length(a.code_pair.logical2_support)
        return collect((n_xstabs + sz_1 + 1):(n_xstabs + sz_1 + w_l + w_r))
    else
        sz_2 = size(parity_matrix_z(a.code_pair.code2), 1)
        w = a.adapter_width
        return collect((n_xstabs + sz_1 + sz_2 + 1):(n_xstabs + sz_1 + sz_2 + 2 * w))
    end
end

"""
Single-logical deformed code (paper Figure 3 / Eq. 9): attach one auxiliary
graph and cycle checks so the Z operator on `logical_support` becomes a
stabilizer product. Use [`build_adapter`](@ref) for joint measurements.

```
H_X_deform = [ H_X | M  ]      original X-stabs deformed by matchings
             [ 0   | N  ]      cycle X-checks (one per basis cycle)

H_Z_deform = [ H_Z | 0   ]     original Z-stabs unchanged
             [ F   | G^T ]     vertex Z-checks (one per port vertex)
```

`max_cycle_len = nothing` defaults to the code's max X-stab weight (paper
§VII.A heuristic). `chords` forwards to [`cellulate_long_cycles!`](@ref).

# Paper Table II reproduction

| Original code | Logical | Deformed params | Default cellulation? |
|---|---|---|:---:|
| BB[[98,6,12]] | Z_1 | [[121,5,12]] | ✓ |
| LP[[200,20,10]] | Z_2 | [[220,19,10]] | ✓ |
| BB[[98,6,12]] (full-rank basis) | Z_3 | [[115,5,12]] | needs `chords=[(1,11)]` |

For BB / Z_3, our default midpoint chord `(4, 9)` and the paper's `(1, 11)`
are both valid 17-edge cellulations but give different single-logical
distance (≤ 11 vs 12). The joint case via [`build_adapter`](@ref) is
unaffected; the bridge X-checks absorb the extra weight-`(d-1)` logical.

```jldoctest
julia> using QuantumClifford, QECCore

julia> using QuantumClifford.ECC

julia> z = logz_ops(Surface(3, 3))[1];

julia> dc = deform_code(Surface(3, 3), z);

julia> (code_n(dc), code_k(dc))
(15, 0)
```
"""
function deform_code(
    code,
    logical;
    max_cycle_len::Union{Int,Nothing} = nothing,
    chords::Union{Nothing, Vector{Tuple{Int,Int}}} = nothing,
)::CSS
    Hx_dense = Matrix{Bool}(parity_matrix_x(code))
    Hz_dense = Matrix{Bool}(parity_matrix_z(code))
    support = _as_support(logical, code, "logical")
    mcl = max_cycle_len === nothing ?
          Int(maximum(sum(Hx_dense; dims = 2))) : max_cycle_len
    aux = cellulate_long_cycles!(
        build_initial_aux_graph(support, code);
        max_cycle_len = mcl, chords = chords)

    n_code = size(Hx_dense, 2)
    sx = size(Hx_dense, 1)
    sz = size(Hz_dense, 1)
    m  = ne(aux.graph)

    M  = stabilizer_modification_matrix(aux, sx)
    N  = aux.cycle_basis_matrix
    F  = port_function_as_matrix(aux, n_code)
    GT = incidence_matrix_transpose(aux)
    n_c = size(N, 1)

    Hx_code = sparse(Hx_dense)
    Hz_code = sparse(Hz_dense)
    Z = (r::Int, c::Int) -> spzeros(Bool, r, c)

    Hx_row1 = hcat(Hx_code,        M)
    Hx_row2 = hcat(Z(n_c, n_code), N)
    Hx_deform = vcat(Hx_row1, Hx_row2)

    Hz_row1 = hcat(Hz_code, Z(sz, m))
    Hz_row2 = hcat(F,       GT)
    Hz_deform = vcat(Hz_row1, Hz_row2)

    CSS(Matrix(Hx_deform), Matrix(Hz_deform))
end
