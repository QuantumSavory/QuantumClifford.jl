using Graphs
using Graphs: Graph, SimpleGraph, add_edge!, neighbors, src, dst, has_edge
using SparseArrays
using LinearAlgebra

"""
Constructs a *Tanner graph* from a given parity-check matrix `H`, where
rows correspond to *check nodes* and columns to *variable nodes*.

The resulting *bipartite* graph indexes variable nodes as `1:n` and check
nodes as `n+1:n+m` for an `m × n` matrix `H`."""
function tanner_graph_from_parity_matrix(H::SparseMatrixCSC{Bool,Int})
    m, n = size(H)
    g = SimpleGraph(n + m)
    var_nodes = collect(1:n)
    check_nodes = collect(n+1:n+m)
    for row in 1:m, col in 1:n
        H[row, col] && add_edge!(g, col, n + row)
    end
    return (graph=g, left=var_nodes, right=check_nodes)
end

"""Generate a random bipartite graph with `n_vars` variable nodes,
`n_checks` check nodes, and edges added with probability `edge_prob`."""
function generate_random_bipartite_graph(n_vars::Int, n_checks::Int, edge_prob::Float64)
    g = SimpleGraph(n_vars + n_checks)
    var_nodes = collect(1:n_vars)
    check_nodes = collect(n_vars+1:n_vars+n_checks)
    for v in var_nodes, c in check_nodes
        rand() < edge_prob && add_edge!(g, v, c)
    end
    return (graph=g, left=var_nodes, right=check_nodes)
end

"""Reconstructs the parity-check matrix from a Tanner graph `g`.

!!! note
    The first block of vertices corresponds to variable nodes and
    the remaining to check nodes. The function returns a parity check
    matrix `H` of size `(number of check nodes) × (number of variable nodes).`
"""
function parity_matrix_from_tanner_graph(g::Graph, var_nodes::Vector{Int}, check_nodes::Vector{Int})
    n_vars = length(var_nodes)
    n_checks = length(check_nodes)
    H = spzeros(Bool, n_checks, n_vars)
    for (check_idx, check_node) in enumerate(check_nodes)
        for neighbor in neighbors(g, check_node)
            neighbor in var_nodes && (H[check_idx, neighbor] = true)
        end
    end
    return H
end

#############################
# Product Tanner Graphs
#############################

function _construct_product_var_nodes(V1, V2, C1, C2)
    var_nodes = Tuple{Symbol,Int,Int}[]
    for v1 in V1, v2 in V2
        push!(var_nodes, (:v, v1, v2))
    end
    for c1 in C1, c2 in C2
        push!(var_nodes, (:c, c1, c2))
    end
    return var_nodes
end

function _construct_product_mappings(var_nodes, check_nodes)
    var_mapping = Dict((node => i for (i, node) in enumerate(var_nodes)))
    check_mapping = Dict((node => i for (i, node) in enumerate(check_nodes)))
    return var_mapping, check_mapping
end

function _add_product_edges!(
    graph, var_nodes, var_mapping, check_mapping,
    G1, G2, V_set, C_set;
    edge_condition, make_check_node
)
    n_vars = length(var_nodes)
    for node in var_nodes
        type, a, b = node
        iter_set = type == :v ? C_set : V_set
        for x in iter_set
            edge_condition(type, G1, G2, a, b, x) || continue
            check_node = make_check_node(type, a, b, x)
            var_idx = var_mapping[node]
            check_idx = check_mapping[check_node]
            add_edge!(graph, var_idx, n_vars + check_idx)
        end
    end
end

"""X-type specific conditions"""
const x_conditions = (
    edge_condition = (type, G1, G2, a, b, x) ->
        type == :v ? has_edge(G1.graph, a, x) : has_edge(G2.graph, x, b),
    make_check_node = (type, a, b, x) ->
        type == :v ? (x, b) : (a, x)
)

"""Z-type specific conditions"""
const z_conditions = (
    edge_condition = (type, G1, G2, a, b, x) -> 
        type == :v ? has_edge(G2.graph, b, x) : has_edge(G1.graph, x, a),
    make_check_node = (type, a, b, x) ->
        type == :v ? (a, x) : (x, b)
)

"""Add edges for X-type Tanner graph.
- For `:v` nodes: iterate over `c` in `C1` and add edge if `G1` has an edge `(v1, c)`.
- For `:c` nodes: iterate over `v` in `V2` and add edge if `G2` has an edge `(v, c2)`."""
add_product_edges_X!(args...) = _add_product_edges!(args...; x_conditions...)

""" Add edges for Z-type Tanner graph.
For `:v` nodes: iterate over `c` in `C2` and add edge if `G2` has an edge `(v2, c)`.
For `:c` nodes: iterate over `v` in `V1` and add edge if `G1` has an edge `(v, c1)`."""
add_product_edges_Z!(args...) = _add_product_edges!(args...; z_conditions...)

"""
Constructs the product Tanner graph for `X-type` checks.

Given Tanner graphs `G₁ = T(V₁, C₁, E₁)` and `G₂ = T(V₂, C₂, E₂)`, define:
- **Variable nodes**: `V = (V₁ × V₂) ∪ (C₁ × C₂)`
- **Check nodes for `X-type`**: `Cₓ = C₁ × V₂`

Edges are added as follows:
- A variable node `(v₁, v₂) ∈ (V₁ × V₂)` connects to a check node
`(c₁, v₂)` if `(v₁, c₁)` is an edge in `G₁`.
- A variable node `(c₁, c₂) ∈ (C₁ × C₂)` connects to a check node
`(c₁, v₂)`if `(v₂, c₂)` is an edge in `G₂`.

The resulting bipartite graph, `G₁ ×ₓ G₂`, defines the classical code `Cₓ`.
"""
function product_tanner_graph_X(G1, G2)
    V1, C1 = G1.left, G1.right
    V2 = G2.left
    var_nodes = _construct_product_var_nodes(V1, V2, C1, G2.right)
    check_nodes = [(c, v) for c in C1 for v in V2]
    var_mapping, check_mapping = _construct_product_mappings(var_nodes, check_nodes)
    total_vertices = length(var_nodes) + length(check_nodes)
    pg = SimpleGraph(total_vertices)
    add_product_edges_X!(pg, var_nodes, var_mapping, check_mapping, G1, G2, V2, C1)
    return (graph=pg, var_nodes=var_nodes, check_nodes=check_nodes,
            var_mapping=var_mapping, check_mapping=check_mapping)
end

"""
Constructs the product Tanner graph for `Z-type` checks.

Given Tanner graphs `G₁ = T(V₁, C₁, E₁)` and `G₂ = T(V₂, C₂, E₂)`, define:
- **Variable nodes**: `V = (V₁ × V₂) ∪ (C₁ × C₂)`
- **Check nodes for `Z-type`**: `C𝓏 = V₁ × C₂`

Edges are added as follows:
- A variable node `(v₁, v₂) ∈ (V₁ × V₂)` connects to a check node
`(v₁, c₂)` if `(v₂, c₂)` is an edge in `G₂`.
- A variable node `(c₁, c₂) ∈ (C₁ × C₂)` connects to a check node
`(v₁, c₂)` if `(v₁, c₁)` is an edge in `G₁`.

The resulting bipartite graph, `G₁ ×𝓏 G₂`, defines the classical code `C𝓏`.
"""
function product_tanner_graph_Z(G1, G2)
    V1, C1 = G1.left, G1.right
    C2 = G2.right
    var_nodes = _construct_product_var_nodes(V1, G2.left, C1, G2.right)
    check_nodes = [(v, c) for v in V1 for c in C2]
    var_mapping, check_mapping = _construct_product_mappings(var_nodes, check_nodes)
    total_vertices = length(var_nodes) + length(check_nodes)
    pg = SimpleGraph(total_vertices)
    add_product_edges_Z!(pg, var_nodes, var_mapping, check_mapping, G1, G2, V1, C2)
    return (graph=pg, var_nodes=var_nodes, check_nodes=check_nodes,
            var_mapping=var_mapping, check_mapping=check_mapping)
end

"""
Extracts the parity-check matrix from a product Tanner graph `PG`.

Vertices ``1\\,..\\,|\\mathit{var\\_nodes}|`` correspond to variable nodes
(matrix columns), while the remaining vertices correspond to check
nodes (matrix rows).

This matrix defines the classical codes `Cₓ = Cₓ(G₁ × G₂)` and `C𝓏 = C𝓏(G₁ × G₂)`,
which together construct the quantum CSS `Q(G₁ × G₂)` code.
"""
function product_parity_matrix_from_tanner_graph(graph::SimpleGraph, var_nodes::Vector, check_nodes::Vector)
    n_vars = length(var_nodes)
    n_checks = length(check_nodes)
    H = spzeros(Bool, n_checks, n_vars)
    for e in edges(graph)
        u, v = src(e), dst(e)
        if u ≤ n_vars && v > n_vars
            H[v - n_vars, u] = true
        elseif v ≤ n_vars && u > n_vars
            H[u - n_vars, v] = true
        end
    end
    return H
end

"""
Constructs a bipartite Tanner graph representing the cycle code
of length `n`.

- The graph consists of `n` variable nodes and `n` check nodes,
totaling `2n` nodes.
- Each check node `i` (for `i` in `1...n`) connects to variable nodes `i`
 and (i+1 \\mod n).

This structure results in a parity-check matrix `H`, where each row contains
exactly two `1`s, encoding the edges of a length-`n` cycle.
"""
function cycle_tanner_graph(n::Int)
    g = SimpleGraph(2n)
    left_nodes = collect(1:n)
    right_nodes = collect(n+1:2n)
    for i in 1:n
        check_node = n + i
        add_edge!(g, i, check_node)
        add_edge!(g, mod1(i+1, n), check_node)
    end
    return (graph=g, left=left_nodes, right=right_nodes)
end

"""
Represents the CSS `Q(G₁ × G₂)` quantum tanner graph product code
constructed from two sparse parity-check matrices H₁ and H₂.

# CSS Code 𝑄(𝐺₁ × 𝐺₂) Associated with a Tanner Graph Product

The CSS code `𝑄(G₁ × G₂)` is constructed from two given sparse parity‐check
matrices `H₁` and `H₂`. This construction follows a specific Tanner graph
product [tillich2013quantum](@cite), defined as follows:

!!! note
    Let `G₁ = 𝑇(V₁, C₁, E₁)` and `G₂ = 𝑇(V₂, C₂, E₂)` be two Tanner graphs. Define the
    vertex and check sets of the product graph as follows: `V = (V₁ × V₂) ∪ (C₁ × C₂)`
    and `C = (C₁ × V₂) ∪ (V₁ × C₂)`. Thus, the product graph `G₁ × G₂` is a bipartite
    graph with vertex set `V ∪ C`.

Two Tanner subgraphs are then defined:

- **G₁ ×ₓ G₂**: A subgraph with variable nodes `V` and check nodes `C₁ × V₂`.
- **G₁ ×𝓏 G₂**: A subgraph with variable nodes `V` and check nodes `V₁ × C₂`.

The union of their edge sets forms `G₁ × G₂`. The classical codes corresponding to
these Tanner graphs are given by: `Cₓ = Cₓ(G₁ × G₂)` and `C𝓏 = C𝓏(G₁ × G₂)`. Together,
these codes define the CSS code 𝑄(G₁ × G₂) quantum tanner graph product code.

# 𝑄(𝐺₁ × 𝐺₂)

The `𝑄(𝐺₁ × 𝐺₂)` quantum LDPC codes represent a broader generalization of **quantum
expander codes** which are derived from the Leverrier-Tillich-Zémor construction [tillich2013quantum](@cite).

```jldoctest examples
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> using SparseArrays; # hide

julia> H1 = sparse(Bool[1 0 1 0; 0 1 0 1; 1 1 0 0]);

julia> H2 = sparse(Bool[1 1 0;0 1 1]);

julia> c = parity_checks(QuantumTannerGraphProduct(H1, H2))
+ X_____X_____X_____
+ _X_____X____XX____
+ __X_____X____X____
+ ___X_____X____X___
+ ____X_____X___XX__
+ _____X_____X___X__
+ X__X____________X_
+ _X__X___________XX
+ __X__X___________X
+ ZZ__________Z___Z_
+ _ZZ__________Z___Z
+ ___ZZ_________Z_Z_
+ ____ZZ_________Z_Z
+ ______ZZ____Z_____
+ _______ZZ____Z____
+ _________ZZ___Z___
+ __________ZZ___Z__
```

# Quantum Expander code

The `𝑄(𝐺₁ × 𝐺₂)` code is more general than the standard quantum expander code
[leverrier2015quantum](@cite) construction. The quantum expander code construction
corresponds to the specific case where `G = G1 = G2`​.

```jldoctest examples
julia> H = sparse(parity_checks(RepCode(3)));

julia> c = parity_checks(QuantumTannerGraphProduct(H, H))
+ X__X_____X_X______
+ _X__X____XX_______
+ __X__X____XX______
+ ___X__X_____X_X___
+ ____X__X____XX____
+ _____X__X____XX___
+ X_____X________X_X
+ _X_____X_______XX_
+ __X_____X_______XX
+ ZZ_______Z_____Z__
+ _ZZ_______Z_____Z_
+ Z_Z________Z_____Z
+ ___ZZ____Z__Z_____
+ ____ZZ____Z__Z____
+ ___Z_Z_____Z__Z___
+ ______ZZ____Z__Z__
+ _______ZZ____Z__Z_
+ ______Z_Z_____Z__Z
```

"""
struct QuantumTannerGraphProduct <: AbstractECC
    H1::AbstractMatrix
    H2::AbstractMatrix
end

function parity_checks_xz(c::QuantumTannerGraphProduct)
    # Build Tanner graphs
    G1 = tanner_graph_from_parity_matrix(c.H1)
    G2 = tanner_graph_from_parity_matrix(c.H2)
    # Construct product graphs
    PG_X = product_tanner_graph_X(G1, G2)
    PG_Z = product_tanner_graph_Z(G1, G2)
    # Extract parity-check matrices
    HX = product_parity_matrix_from_tanner_graph(PG_X.graph, PG_X.var_nodes, PG_X.check_nodes)
    HZ = product_parity_matrix_from_tanner_graph(PG_Z.graph, PG_Z.var_nodes, PG_Z.check_nodes)
    return Matrix{Bool}(HX), Matrix{Bool}(HZ)
end

""" Constructs a 𝑄(𝐺₁ × 𝐺₂) quantum Tanner graph product code
using cyclic Tanner graphs of length `2m`.

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> m = 10;

julia> c = parity_checks(CyclicQuantumTannerGraphProduct(m));

julia> code_n(c), code_k(c)
(800, 2)
```
"""
struct CyclicQuantumTannerGraphProduct <: AbstractECC
    m::Int
    function CyclicQuantumTannerGraphProduct(m::Int)
        m > 0 || throw(ArgumentError("m must be positive"))
        new(m)
    end
end

function parity_checks_xz(Q::CyclicQuantumTannerGraphProduct)
    n = 2 * Q.m
    G1 = cycle_tanner_graph(n)
    G2 = cycle_tanner_graph(n)
    PG_X = product_tanner_graph_X(G1, G2)
    PG_Z = product_tanner_graph_Z(G1, G2)
    HX = product_parity_matrix_from_tanner_graph(PG_X.graph, PG_X.var_nodes, PG_X.check_nodes)
    HZ = product_parity_matrix_from_tanner_graph(PG_Z.graph, PG_Z.var_nodes, PG_Z.check_nodes)
    return (Matrix{Bool}(HX), Matrix{Bool}(HZ))
end

iscss(::Type{<:Union{QuantumTannerGraphProduct,CyclicQuantumTannerGraphProduct}}) = true

parity_checks_x(c::Union{QuantumTannerGraphProduct,CyclicQuantumTannerGraphProduct}) = parity_checks_xz(c)[1]
parity_checks_z(c::Union{QuantumTannerGraphProduct,CyclicQuantumTannerGraphProduct}) = parity_checks_xz(c)[2]

function parity_checks(c::Union{QuantumTannerGraphProduct,CyclicQuantumTannerGraphProduct})
    hx, hz = parity_checks_xz(c)
    return parity_checks(CSS(parity_checks_xz(c)...))
end
