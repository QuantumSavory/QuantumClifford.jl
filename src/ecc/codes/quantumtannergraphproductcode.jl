using Graphs
using SparseArrays
using LinearAlgebra

"""
Represents a bipartite Tanner graph used in classical coding theory. It
connects *variable nodes* (codeword bits or qubits) to *check nodes*
(parity constraints) and serves as the foundation for constructing
product graphs in quantum CSS codes.
"""
struct BipartiteGraph
    """The underlying undirected graph (from `Graphs.jl`)."""
    g::Graph
    """Indices of variable nodes."""
    left_nodes::Vector{Int}
    """Indices of check nodes."""
    right_nodes::Vector{Int}
end

"""
Constructs a Tanner graph from a given parity-check matrix `H`, where
rows correspond to check nodes and columns to variable nodes. The resulting
`BipartiteGraph` indexes variable nodes as `1:n` and check nodes as `n+1:n+m`
for an `m × n` matrix `H`. 

!!! note The input `H` must be a sparse Boolean matrix representing the
parity-check matrix.
"""
function tanner_graph_from_parity_matrix(H::SparseMatrixCSC{Bool,Int})
    m, n = size(H)
    g = Graph(n + m)  # Total nodes: n variable nodes + m check nodes
    left_nodes = collect(1:n)
    right_nodes = collect(n+1:n+m)
    # For every 1 in H, add an edge between the corresponding variable and check nodes.
    for row in 1:m
        for col in 1:n
            if H[row, col]
                add_edge!(g, col, n + row)
            end
        end
    end
    return BipartiteGraph(g, left_nodes, right_nodes)
end

"""
Reconstructs the parity-check matrix from a Tanner graph `bg`, assuming
the first block of vertices corresponds to variable nodes and the remaining
to check nodes. The function returns a sparse Boolean matrix `H` of size
`(number of check nodes) × (number of variable nodes).`

!!! note The input `bg` must be a `BipartiteGraph` instance.
"""
function parity_matrix_from_tanner_graph(bg::BipartiteGraph)
    n_vars = length(bg.left_nodes)
    n_checks = length(bg.right_nodes)
    H = spzeros(Bool, n_checks, n_vars)
    # For each check node, mark a 1 in H for each incident variable node.
    for check_node in bg.right_nodes
        check_idx = check_node - n_vars
        for neighbor in neighbors(bg.g, check_node)
            if neighbor in bg.left_nodes
                var_idx = neighbor
                H[check_idx, var_idx] = true
            end
        end
    end
    return H
end

#############################
# Product Tanner Graphs
#############################

"""
Represents a product Tanner graph formed by combining two Tanner graphs, `G₁`
and `G₂`, to define `X`-type or `Z`-type checks. The variable node set is:
`V = (V₁ × V₂) ∪ (C₁ × C₂)`. The check node set depends on the type of product
graph:
- **X-type**: `Cₓ = C₁ × V₂`
- **Z-type**: `C𝓏 = V₁ × C₂`
"""
struct ProductTannerGraph
    """The underlying graph, where vertices `1..|var_nodes|` correspond to variable nodes,
       and `|var_nodes|+1..end` correspond to check nodes."""
    g::Graph
    """A list of variable nodes represented as `(type, index1, index2)`."""
    var_nodes::Vector{Tuple{Symbol,Int,Int}}
    """A list of check nodes represented as `(index1, index2)`."""
    check_nodes::Vector{Tuple{Int,Int}}
    """Maps a variable node label `(type, index1, index2)` to its corresponding vertex index."""
    var_mapping::Dict{Tuple{Symbol,Int,Int},Int}
    """Maps a check node label `(index1, index2)` to its corresponding vertex index."""
    check_mapping::Dict{Tuple{Int,Int},Int}
end

"""
Constructs the product Tanner graph for `X-type` checks.

Given Tanner graphs `G₁ = T(V₁, C₁, E₁)` and `G₂ = T(V₂, C₂, E₂)`, define:

- **Variable nodes**: `V = (V₁ × V₂) ∪ (C₁ × C₂)`
- **Check nodes for X-type**: `Cₓ = C₁ × V₂`

Edges are added as follows:
1. A variable node `(v₁, v₂) ∈ (V₁ × V₂)` connects to a check node `(c₁, v₂)`
   if `(v₁, c₁)` is an edge in `G₁`.
2. A variable node `(c₁, c₂) ∈ (C₁ × C₂)` connects to a check node `(c₁, v₂)`
   if `(v₂, c₂)` is an edge in `G₂`.

The resulting bipartite graph, `G₁ ×ₓ G₂`, defines the classical code `Cₓ`.
"""
function product_tanner_graph_X(G1::BipartiteGraph, G2::BipartiteGraph)
    V1 = G1.left_nodes
    C1 = G1.right_nodes
    V2 = G2.left_nodes
    C2 = G2.right_nodes
    # Construct variable node list: (V₁ × V₂) ∪ (C₁ × C₂).
    var_nodes = Tuple{Symbol,Int,Int}[]
    for v1 in V1, v2 in V2
        push!(var_nodes, (:v, v1, v2))
    end
    for c1 in C1, c2 in C2
        push!(var_nodes, (:c, c1, c2))
    end
    # Define check node list for X-type: Cₓ = C₁ × V₂.
    check_nodes = [(c1, v2) for c1 in C1 for v2 in V2]
    # Create mappings.
    var_mapping = Dict{Tuple{Symbol,Int,Int},Int}()
    for (i, node) in enumerate(var_nodes)
        var_mapping[node] = i
    end
    check_mapping = Dict{Tuple{Int,Int},Int}()
    for (i, node) in enumerate(check_nodes)
        check_mapping[node] = i
    end
    total_vertices = length(var_nodes) + length(check_nodes)
    pg = Graph(total_vertices)
    # Add edges for X-type checks 
    # For variable nodes of type (:v, v1, v2).
    for node in var_nodes
        if node[1] == :v
            (_, v1, v2) = node
            for c1 in C1
                check_node = (c1, v2)
                if has_edge(G1.g, v1, c1)
                    var_index = var_mapping[node]
                    check_index = check_mapping[check_node]
                    add_edge!(pg, var_index, length(var_nodes) + check_index)
                end
            end
        end
    end
    # For variable nodes of type (:c, c1, c2).
    for node in var_nodes
        if node[1] == :c
            (_, c1, c2) = node
            for v2 in V2
                check_node = (c1, v2)
                if has_edge(G2.g, v2, c2)
                    var_index = var_mapping[node]
                    check_index = check_mapping[check_node]
                    add_edge!(pg, var_index, length(var_nodes) + check_index)
                end
            end
        end
    end
    return ProductTannerGraph(pg, var_nodes, check_nodes, var_mapping, check_mapping)
end

"""
Constructs the product Tanner graph for `Z-type` checks.

Given Tanner graphs `G₁ = T(V₁, C₁, E₁)` and `G₂ = T(V₂, C₂, E₂)`, define:

- **Variable nodes**: `V = (V₁ × V₂) ∪ (C₁ × C₂)`
- **Check nodes for `Z-type`**: `C𝓏 = V₁ × C₂`

Edges are added as follows:
1. Connect a variable node `(v₁, v₂) ∈ (V₁ × V₂)` to a check node `(v₁, c₂)`
   if `(v₂, c₂)` is an edge in `G₂`
2. Connect a variable node `(c₁, c₂) ∈ (C₁ × C₂)` to a check node `(v₁, c₂)`
   if `(v₁, c₁)` is an edge in `G₁`.

The resulting bipartite graph, `G₁ ×𝓏 G₂`, defines the classical code `C𝓏`.
"""
function product_tanner_graph_Z(G1::BipartiteGraph, G2::BipartiteGraph)
    V1 = G1.left_nodes
    C1 = G1.right_nodes
    V2 = G2.left_nodes
    C2 = G2.right_nodes
    # Construct variable node list: (V₁ × V₂) ∪ (C₁ × C₂).
    var_nodes = Tuple{Symbol,Int,Int}[]
    for v1 in V1, v2 in V2
        push!(var_nodes, (:v, v1, v2))
    end
    for c1 in C1, c2 in C2
        push!(var_nodes, (:c, c1, c2))
    end
    # Define check node list for Z-type: C𝓏 = V₁ × C₂.
    check_nodes = [(v1, c2) for v1 in V1 for c2 in C2]
    # Create mappings.
    var_mapping = Dict{Tuple{Symbol,Int,Int},Int}()
    for (i, node) in enumerate(var_nodes)
        var_mapping[node] = i
    end
    check_mapping = Dict{Tuple{Int,Int},Int}()
    for (i, node) in enumerate(check_nodes)
        check_mapping[node] = i
    end
    total_vertices = length(var_nodes) + length(check_nodes)
    pg = Graph(total_vertices)
    # Add edges for Z-type checks
    # For variable nodes of type (:v, v1, v2).
    for node in var_nodes
        if node[1] == :v
            (_, v1, v2) = node
            for c2 in C2
                check_node = (v1, c2)
                if has_edge(G2.g, v2, c2)
                    var_index = var_mapping[node]
                    check_index = check_mapping[check_node]
                    add_edge!(pg, var_index, length(var_nodes) + check_index)
                end
            end
        end
    end
    # For variable nodes of type (:c, c1, c2).
    for node in var_nodes
        if node[1] == :c
            (_, c1, c2) = node
            for v1 in V1
                check_node = (v1, c2)
                if has_edge(G1.g, v1, c1)
                    var_index = var_mapping[node]
                    check_index = check_mapping[check_node]
                    add_edge!(pg, var_index, length(var_nodes) + check_index)
                end
            end
        end
    end
    return ProductTannerGraph(pg, var_nodes, check_nodes, var_mapping, check_mapping)
end

"""
Extracts the parity-check matrix from a product Tanner graph `PG`.

Vertices `1..|var_nodes|` correspond to variable nodes (matrix columns),
while the remaining vertices correspond to check nodes (matrix rows).

This matrix defines the classical codes `Cₓ = Cₓ(G₁ × G₂)` and `C𝓏 = C𝓏(G₁ × G₂)`,
which together construct the CSS code `Q(G₁ × G₂)`.
"""
function product_parity_matrix_from_tanner_graph(PG::ProductTannerGraph)
    n_vars = length(PG.var_nodes)
    n_checks = length(PG.check_nodes)
    H = spzeros(Bool, n_checks, n_vars)
    for e in edges(PG.g)
        u, v = src(e), dst(e)
        if u ≤ n_vars && v > n_vars
            row = v - n_vars
            col = u
            H[row, col] = true
        elseif v ≤ n_vars && u > n_vars
            row = u - n_vars
            col = v
            H[row, col] = true
        end
    end
    return H
end

"""
Checks whether the parity-check matrices H_X and H_Z are orthogonal in GF(2),
a key requirement for CSS code construction. In product Tanner graphs, the
existence of 4-cycles ensures orthogonality. 
"""
function verify_orthogonality(HX::SparseMatrixCSC{Bool,Int}, HZ::SparseMatrixCSC{Bool,Int})
    product_mod2 = mod.(Int.(HX) * transpose(Int.(HZ)), 2)
    return all(product_mod2 .== 0)
end

"""
Constructs a bipartite Tanner graph representing the cycle code
of length `n`.

- The graph consists of `n` variable nodes and `n` check nodes,
totaling `2n` nodes.
- Each check node i (for i in 1..n) connects to variable nodes i
 and (i+1 \\mod n).

This structure results in a parity-check matrix  H, where each row contains
exactly two `1`s, encoding the edges of a length-n cycle.
"""
function cycle_tanner_graph(n::Int)
    # Total nodes = n variables + n checks.
    g = Graph(2n)
    left_nodes  = collect(1:n)
    right_nodes = collect(n+1:2n)
    for i in 1:n
        # The i-th check node is node index n+i
        check_node = n + i
        # Connect it to variable i and variable (i mod n) + 1
        add_edge!(g, i, check_node)
        add_edge!(g, ((i % n) + 1), check_node)
    end
    return BipartiteGraph(g, left_nodes, right_nodes)
end

"""
Represents the CSS code Q(G₁ × G₂) constructed from two sparse parity-check
matrices H₁ and H₂.

# CSS Code 𝑄(𝐺₁ × 𝐺₂) Associated with a Graph Product

The CSS code 𝑄(G₁ × G₂) is constructed from two given sparse parity‐check
matrices H₁ and H₂. This construction follows a specific Tanner graph
product, defined as follows:

Let G₁ = 𝑇(V₁, C₁, E₁) and G₂ = 𝑇(V₂, C₂, E₂) be two Tanner graphs. Define the
vertex and check sets of the product graph as follows: V = (V₁ × V₂) ∪ (C₁ × C₂)
and C = (C₁ × V₂) ∪ (V₁ × C₂) Thus, the product graph G₁ × G₂ is a bipartite graph
with vertex set V ∪ C.

Two Tanner subgraphs are then defined:

- **G₁ ×ₓ G₂**: A subgraph with variable nodes V and check nodes C₁ × V₂.
- **G₁ ×𝓏 G₂**: A subgraph with variable nodes V and check nodes V₁ × C₂.

The union of their edge sets forms G₁ × G₂. The classical codes corresponding to
these Tanner graphs are given by: Cₓ = Cₓ(G₁ × G₂) and C𝓏 = C𝓏(G₁ × G₂). Together,
these codes define the CSS code 𝑄(G₁ × G₂).

## Construction Method

Given two parity-check matrices H₁ and H₂, the CSS code 𝑄(G₁ × G₂) is constructed
through the following steps:

- Convert H₁ and H₂ into Tanner graphs G₁ and G₂ using `tanner_graph_from_parity_matrix`.
- Construct the product Tanner graphs for the X-type and Z-type checks using
`product_tanner_graph_X` and `product_tanner_graph_Z`.
- Extract the final parity‐check matrices Hₓ and H𝓏 from these graphs via
`product_parity_matrix_from_tanner_graph`.

Since the Tanner graphs can be fully reconstructed from H₁ and H₂, only these input
matrices and the resulting parity-check matrices Hₓ and H𝓏 need to be stored in order
to define the CSS code.

# Example

```jldoctest
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
"""
struct QuantumTannerGraphProduct <: AbstractECC
    H1::SparseMatrixCSC{Bool,Int}
    H2::SparseMatrixCSC{Bool,Int}

    function QuantumTannerGraphProduct(H1::SparseMatrixCSC, H2::SparseMatrixCSC)
        isa(H1, SparseMatrixCSC{Bool,Int}) || throw(ArgumentError("H1 must be a sparse Boolean matrix."))
        isa(H2, SparseMatrixCSC{Bool,Int}) || throw(ArgumentError("H2 must be a sparse Boolean matrix."))
        new(H1, H2)
    end
end

function iscss(::Type{QuantumTannerGraphProduct})
    return true
end

function parity_checks_xz(c::QuantumTannerGraphProduct)
    # Build Tanner graphs from H1 and H2.
    G1 = tanner_graph_from_parity_matrix(c.H1)
    G2 = tanner_graph_from_parity_matrix(c.H2)
    # Construct product Tanner graphs for X-type and Z-type checks.
    PG_X = product_tanner_graph_X(G1, G2)
    PG_Z = product_tanner_graph_Z(G1, G2)
    # Extract the final parity-check matrices.
    HX = product_parity_matrix_from_tanner_graph(PG_X)
    HZ = product_parity_matrix_from_tanner_graph(PG_Z)
    hx = Matrix{Bool}(HX)
    hz = Matrix{Bool}(HZ)
    return hx, hz
end

parity_checks_x(c::QuantumTannerGraphProduct) = parity_checks_xz(c)[1]

parity_checks_z(c::QuantumTannerGraphProduct) = parity_checks_xz(c)[2]

parity_checks(c::QuantumTannerGraphProduct) = parity_checks(CSS(parity_checks_xz(c)...))

""" Constructs a 𝑄(𝐺₁ × 𝐺₂) quantum Tanner graph product code using cyclic
Tanner graphs of length 2m.

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
        m > 0 || throw(ArgumentError("m must be a positive integer."))
        return new(m)
    end
end

function iscss(::Type{CyclicQuantumTannerGraphProduct})
    return true
end

function parity_checks_xz(Q::CyclicQuantumTannerGraphProduct)
    n = 2 * Q.m
    G1 = cycle_tanner_graph(n)
    G2 = cycle_tanner_graph(n)
    PG_X = product_tanner_graph_X(G1, G2)
    PG_Z = product_tanner_graph_Z(G1, G2)
    HX = product_parity_matrix_from_tanner_graph(PG_X)
    HZ = product_parity_matrix_from_tanner_graph(PG_Z)
    hx = Matrix{Bool}(HX)
    hz = Matrix{Bool}(HZ)
    return hx, hz
end

parity_checks_x(c::CyclicQuantumTannerGraphProduct) = parity_checks_xz(c)[1]

parity_checks_z(c::CyclicQuantumTannerGraphProduct) = parity_checks_xz(c)[2]

parity_checks(c::CyclicQuantumTannerGraphProduct) = parity_checks(CSS(parity_checks_xz(c)...))
