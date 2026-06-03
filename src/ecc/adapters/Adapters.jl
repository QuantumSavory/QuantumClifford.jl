"""
    QuantumClifford.ECC.Adapters

Universal-adapter construction between CSS qLDPC codeblocks via the SkipTree
algorithm of [swaroop2026universal](@cite). See [`build_adapter`](@ref).
"""
module Adapters

using SparseArrays: SparseMatrixCSC, sparse, nnz, nonzeros, rowvals, nzrange,
                    findnz, spzeros
using LinearAlgebra: LinearAlgebra
using Graphs: Graphs, SimpleGraph, AbstractGraph, add_edge!, has_edge, neighbors,
              kruskal_mst, nv, ne, edges, is_connected, connected_components,
              cycle_basis, src, dst, vertices

using DocStringExtensions: TYPEDEF, TYPEDFIELDS

using QECCore: CSS, AbstractCSSCode
import QECCore: code_n, code_s, parity_matrix_x, parity_matrix_z

include("types.jl")
include("skiptree.jl")
include("aux_graph.jl")
include("merge.jl")

export
    build_adapter, build_adapter_intercode, build_adapter_intracode, deform_code,
    skiptree,
    CodePair, AuxiliaryGraph, SkipTreeOutput, Adapter, AdapterMergedCode

end # module Adapters
