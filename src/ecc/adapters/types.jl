"""
$TYPEDEF

A pair of logical Z operators we want to jointly measure via a universal adapter.
`code1` and `code2` are CSS code blocks; `code1 === code2` triggers the intra-code
path. `logical{1,2}_support` are sorted 1-indexed lists of qubits in their respective
codes; each must have weight ≥ d (code distance).

X-type and Y-type joint measurements reduce to the Z-type case via local Cliffords
(paper §IV), so only Z is supported directly.

### Fields
$TYPEDFIELDS
"""
struct CodePair
    "Left CSS codeblock."
    code1::CSS
    "Right CSS codeblock (`=== code1` for intra-code measurement)."
    code2::CSS
    "Sorted 1-indexed qubit indices supporting Z̄_l in code1."
    logical1_support::Vector{Int}
    "Sorted 1-indexed qubit indices supporting Z̄_r in code2."
    logical2_support::Vector{Int}
end

"""
$TYPEDEF

Auxiliary graph G for measuring a single logical, built by
[`build_initial_aux_graph`](@ref) and finalised by
[`cellulate_long_cycles!`](@ref). One graph per side of the adapter.

`port_function` is set-valued per Appendix B of [swaroop2026universal](@cite) to
support shared-qubit intra-code measurements; for the usual injective case each
inner vector has length 1.

`stabilizer_matchings[s]` holds vertex pairs `(min, max)` for the perfect matching
μ(L_s) rather than edge indices: `SimpleGraph`'s edge iteration shifts when
cellulation adds chords, but vertex pairs are stable. Translate to current column
index via [`edge_pair_to_index`](@ref).

### Fields
$TYPEDFIELDS
"""
struct AuxiliaryGraph
    "The auxiliary graph G = (V, E)."
    graph::SimpleGraph{Int}
    "`port_function[q]` lists the vertices attached to logical-support qubit `q`."
    port_function::Vector{Vector{Int}}
    "Global qubit index of the ancilla edge-qubit for each edge."
    edge_to_qubit::Vector{Int}
    "Cycle basis matrix N (n_cycles × n_edges); 1 iff edge in basis cycle. Empty until cellulation runs."
    cycle_basis_matrix::SparseMatrixCSC{Bool, Int}
    "Sorted 1-indexed qubits in the codeblock supporting the target logical."
    logical_support::Vector{Int}
    "Vertex pairs `(min, max)` in μ(L_s) per X-stabilizer row; empty rows for non-overlapping stabs."
    stabilizer_matchings::Vector{Vector{Tuple{Int,Int}}}
end

"""
$TYPEDEF

Output of the SkipTree algorithm (Algorithm 2, Appendix E of
[swaroop2026universal](@cite)). Satisfies `T · G · P ≡ H_R(n) (mod 2)` where
`G` is the `m × n` edge-vertex incidence of the graph and `H_R(n)` is the
banded `(n-1) × n` repetition check matrix.

`T` is `(3, 2)`-sparse: row weight ≤ 3, column weight ≤ 2.

### Fields
$TYPEDFIELDS
"""
struct SkipTreeOutput
    "Sparse `(n-1, m)` basis-change matrix; `(3, 2)`-sparse over GF(2)."
    T::SparseMatrixCSC{Bool, Int}
    "Sparse `(n, n)` vertex permutation."
    P::SparseMatrixCSC{Bool, Int}
    "Algorithm target. Only `:H_R` is supported (Algorithm 2)."
    target::Symbol
end

"""
$TYPEDEF

Universal adapter between two CSS codes — the full output of [`build_adapter`](@ref).
`merged_code` is the resulting deformed CSS code in which Z̄_l Z̄_r becomes a product
of stabilizers; the remaining fields preserve the construction for downstream consumers
(decoders, hardware-aware placement, distance specialisations).

`adapter_width = min(|logical1_support|, |logical2_support|)` is the bridge size
per Lemma 9.

### Fields
$TYPEDFIELDS
"""
struct Adapter
    "Input pair of logicals."
    code_pair::CodePair
    "Left auxiliary graph (for `logical1_support`)."
    aux_l::AuxiliaryGraph
    "Right auxiliary graph (for `logical2_support`)."
    aux_r::AuxiliaryGraph
    "SkipTree output for `aux_l`."
    skiptree_l::SkipTreeOutput
    "SkipTree output for `aux_r`."
    skiptree_r::SkipTreeOutput
    "Bridge width |A| = min(|L_l|, |L_r|)."
    adapter_width::Int
    "The merged CSS code containing Z̄_l Z̄_r as a stabilizer product."
    merged_code::CSS
end

"""
$TYPEDEF

`AbstractCSSCode` wrapper around an [`Adapter`](@ref) so the merged code plugs into
the standard `QECCore.AbstractECC` dispatch (`code_n`, `parity_matrix_x`, …) while
still exposing the originating `Adapter` via `.adapter` for consumers that need the
auxiliary graphs or SkipTree outputs.

# Example

```jldoctest
julia> using QuantumClifford, QECCore

julia> using QuantumClifford.ECC: Surface, logz_ops, code_n, code_k, code_s

julia> using QuantumClifford.ECC: parity_matrix_x, parity_matrix_z

julia> using QuantumClifford: stab_to_gf2

julia> c1 = CSS(Matrix{Bool}(parity_matrix_x(Surface(3,3))), Matrix{Bool}(parity_matrix_z(Surface(3,3))));

julia> c2 = CSS(Matrix{Bool}(parity_matrix_x(Surface(3,3))), Matrix{Bool}(parity_matrix_z(Surface(3,3))));

julia> z = sort(findall(!iszero, stab_to_gf2(logz_ops(Surface(3,3)))[1, 14:26]));

julia> wrap = AdapterMergedCode(build_adapter(CodePair(c1, c2, z, z)));

julia> (code_n(wrap), code_k(wrap), code_s(wrap))
(33, 1, 32)
```

### Fields
$TYPEDFIELDS
"""
struct AdapterMergedCode <: AbstractCSSCode
    "The originating universal-adapter assembly."
    adapter::Adapter
end
