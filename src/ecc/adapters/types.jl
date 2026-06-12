"""
Input to [`build_adapter`](@ref): two CSS code blocks plus the qubit supports
of the two logical Z operators to be jointly measured. `code1 === code2`
selects the intra-code path. X- and Y-type joint measurements reduce to the
Z-type case via local Cliffords (paper §IV).

- `code1`, `code2`: the input codes
- `logical1_support`, `logical2_support`: qubit indices supporting Z̄_l and Z̄_r
"""
struct CodePair
    code1::CSS
    code2::CSS
    logical1_support::Vector{Int}
    logical2_support::Vector{Int}
end

"""
Auxiliary graph G constructed for measuring a single logical, built by
[`build_initial_aux_graph`](@ref) and finalised by [`cellulate_long_cycles!`](@ref).
One per side of the adapter.

`port_function[i]` is the list of graph vertices that the i-th qubit of
`logical_support` maps to. In the worked paper examples this list always has
length 1 (injective port function: one qubit → one vertex); the nested-vector
type exists because paper Appendix B generalises to allow one support qubit
to map to several vertices.

`stabilizer_matchings` holds vertex pairs `(min, max)` rather than edge indices,
because `SimpleGraph`'s edge iteration order shifts when cellulation adds chords
and vertex pairs are stable under that.
"""
struct AuxiliaryGraph
    graph::SimpleGraph{Int}
    port_function::Vector{Vector{Int}}
    edge_to_qubit::Vector{Int}
    cycle_basis_matrix::SparseMatrixCSC{Bool, Int}
    logical_support::Vector{Int}
    stabilizer_matchings::Vector{Vector{Tuple{Int,Int}}}
end

"""
Output of the SkipTree algorithm (paper Algorithm 2, Appendix E). Holds the
`(n-1, m)` basis-change matrix and the `n × n` vertex permutation such that
`T · G · P ≡ H_R(n) (mod 2)` with `G` the `m × n` edge-vertex incidence and
`H_R(n)` the banded `(n-1) × n` repetition check matrix. `T` is `(3, 2)`-sparse.
Only the `:H_R` target is supported.
"""
struct SkipTreeOutput
    T::SparseMatrixCSC{Bool, Int}
    P::SparseMatrixCSC{Bool, Int}
    target::Symbol
end

"""
A universal adapter between two CSS codes — the result of [`build_adapter`](@ref).
Subtypes `AbstractCSSCode`, so it plugs directly into the standard `QECCore`
dispatch (`code_n`, `parity_matrix_x`, distance, decoder harnesses) without an
intermediate wrapper.

# Example

```jldoctest
julia> using QuantumClifford, QECCore

julia> using QuantumClifford.ECC: Surface, logz_ops, code_n, code_k, code_s

julia> using QuantumClifford.ECC: parity_matrix_x, parity_matrix_z

julia> using QuantumClifford: stab_to_gf2

julia> c1 = CSS(Matrix{Bool}(parity_matrix_x(Surface(3,3))), Matrix{Bool}(parity_matrix_z(Surface(3,3))));

julia> c2 = CSS(Matrix{Bool}(parity_matrix_x(Surface(3,3))), Matrix{Bool}(parity_matrix_z(Surface(3,3))));

julia> z = sort(findall(!iszero, stab_to_gf2(logz_ops(Surface(3,3)))[1, 14:26]));

julia> adapter = build_adapter(CodePair(c1, c2, z, z));

julia> (code_n(adapter), code_k(adapter), code_s(adapter))
(33, 1, 32)
```
"""
struct Adapter <: AbstractCSSCode
    code_pair::CodePair
    aux_l::AuxiliaryGraph
    aux_r::AuxiliaryGraph
    skiptree_l::SkipTreeOutput
    skiptree_r::SkipTreeOutput
    adapter_width::Int
    merged_code::CSS
end
