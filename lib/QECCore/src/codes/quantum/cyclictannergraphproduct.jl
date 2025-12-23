"""
Represents a *bipartite* Tanner graph used in classical coding theory. It connects *variable
nodes* (codeword bits or qubits) to *check nodes* (parity constraints) and serves as the
foundation for constructing Tanner graphs in quantum CSS codes.
"""
struct BipartiteGraph
    """The underlying undirected graph (from `Graphs.jl`)."""
    g::Graph
    """Indices of variable nodes."""
    var_nodes::Vector{Int}
    """Indices of check nodes."""
    check_nodes::Vector{Int}
end

"""
Constructs a *Tanner graph* from a given parity-check matrix `H`, where rows correspond to *check nodes*
and columns to *variable nodes*. The resulting *bipartite* graph indexes variable nodes as `1:n` and check
nodes as `n+1:n+m` for an `m × n` matrix `H`."""
function tanner_graph_from_parity_matrix(H::AbstractMatrix)
    m, n = size(H)
    g = SimpleGraph(n + m)
    var_nodes = collect(1:n)
    check_nodes = collect(n+1:n+m)
    for row in 1:m, col in 1:n
        H[row, col] && add_edge!(g, col, n + row)
    end
    return BipartiteGraph(g, var_nodes, check_nodes)
end

"""Reconstructs the parity-check matrix from a Tanner graph `g`.

!!! note
    The first block of vertices corresponds to variable nodes and the remaining to check nodes. It
    returns a parity check matrix `H` of size `(number of check nodes) × (number of variable nodes).`
"""
function parity_matrix_from_tanner_graph(bg::BipartiteGraph)
    n_vars = length(bg.var_nodes)
    n_checks = length(bg.check_nodes)
    H = spzeros(Bool, n_checks, n_vars)
    for (check_idx, check_node) in enumerate(bg.check_nodes)
        for neighbor in neighbors(bg.g, check_node)
            neighbor in bg.var_nodes && (H[check_idx, neighbor] = true)
        end
    end
    return H
end

"""
Constructs a bipartite Tanner graph representing the cycle code of length `n`. The graph
consists of `n` variable nodes and `n` check nodes, totaling `2n` nodes. Each check node `i`
(for `i` in `1...n`) connects to variable nodes `i` and (``i+1 \\mod n``). This results in a
parity-check matrix `H`, where each row contains exactly two `1`s, encoding the edges of a
length-`n` cycle.
"""
function cycle_tanner_graph(n::Int)
    g = SimpleGraph(2n)
    var_nodes = collect(1:n)
    check_nodes = collect(n+1:2n)
    for i in 1:n
        check_node = n + i
        add_edge!(g, i, check_node)
        add_edge!(g, mod1(i+1, n), check_node)
    end
    return BipartiteGraph(g, var_nodes, check_nodes)
end

"""
    $TYPEDEF

Constructs a `𝑄(𝐺₁ × 𝐺₂)` quantum Tanner graph product code using cyclic Tanner graphs of length `2m`.

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC;

julia> m = 1;

julia> c = parity_checks(CyclicQuantumTannerGraphProduct(m))
+ X_X_XX__
+ _X_XXX__
+ X_X___XX
+ _X_X__XX
+ ZZ__Z_Z_
+ ZZ___Z_Z
+ __ZZZ_Z_
+ __ZZ_Z_Z

julia> code_n(c), code_k(c)
(8, 2)

julia> m = 10;

julia> c = parity_checks(CyclicQuantumTannerGraphProduct(m));

julia> code_n(c), code_k(c)
(800, 2)
```

### Fields
    $TYPEDFIELDS
"""
struct CyclicQuantumTannerGraphProduct <: AbstractCSSCode
    m::Int
    function CyclicQuantumTannerGraphProduct(m::Int)
        m > 0 || throw(ArgumentError("m must be positive."))
        new(m)
    end
end

function parity_matrix_xz(Q::CyclicQuantumTannerGraphProduct)
    n = 2*Q.m
    G1 = cycle_tanner_graph(n)
    G2 = cycle_tanner_graph(n)
    H1 = parity_matrix_from_tanner_graph(G1)
    H2 = parity_matrix_from_tanner_graph(G2)
    return hgp(H1, H2)
end

parity_matrix_x(c::CyclicQuantumTannerGraphProduct) = parity_matrix_xz(c)[1]

parity_matrix_z(c::CyclicQuantumTannerGraphProduct) = parity_matrix_xz(c)[2]
