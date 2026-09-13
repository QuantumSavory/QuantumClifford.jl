"""
Input bundle for [`build_adapter`](@ref): two CSS codeblocks plus the qubit
supports of the two ``Z`` logicals being jointly measured. The `intracode`
flag indicates whether the two logicals live on the same codeblock (and the
merged code has only one copy of the data qubits) or on two independent
codeblocks. ``X``- and ``Y``-type joint measurements reduce to the ``Z`` case
via local Cliffords (paper §IV).

`code1` and `code2` only need to support `parity_matrix_x` and `parity_matrix_z`
— the constructor verifies both at build time, so any `AbstractECC` (or
duck-typed code) that implements those methods works.
"""
struct CodePair{C1, C2}
    code1::C1
    code2::C2
    logical1_support::Vector{Int}
    logical2_support::Vector{Int}
    intracode::Bool

    function CodePair(code1::C1, code2::C2,
                       logical1_support::AbstractVector{<:Integer},
                       logical2_support::AbstractVector{<:Integer};
                       intracode::Bool = false) where {C1, C2}
        _verify_css_methods(code1, "code1")
        _verify_css_methods(code2, "code2")
        new{C1, C2}(code1, code2,
                    collect(Int, logical1_support),
                    collect(Int, logical2_support),
                    intracode)
    end
end

function _verify_css_methods(code, label::String)
    applicable(parity_matrix_x, code) || throw(ArgumentError(
        "CodePair: $label of type $(typeof(code)) does not implement parity_matrix_x"))
    applicable(parity_matrix_z, code) || throw(ArgumentError(
        "CodePair: $label of type $(typeof(code)) does not implement parity_matrix_z"))
    return nothing
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
Subtypes `AbstractCSSCode`, so it behaves like and can be used as any other
stabilizer code in this library.

# Example

```jldoctest
julia> using QuantumClifford, QECCore

julia> using QuantumClifford.ECC

julia> c1 = Surface(3, 3);

julia> c2 = Surface(3, 3);

julia> z = logz_ops(Surface(3, 3))[1];

julia> adapter = build_adapter(c1, c2, z, z);

julia> (code_n(adapter), code_k(adapter), code_s(adapter))
(33, 1, 32)
```
"""
struct Adapter{C1, C2} <: AbstractCSSCode
    code_pair::CodePair{C1, C2}
    aux_l::AuxiliaryGraph
    aux_r::AuxiliaryGraph
    skiptree_l::SkipTreeOutput
    skiptree_r::SkipTreeOutput
    adapter_width::Int
    merged_code::CSS
end
