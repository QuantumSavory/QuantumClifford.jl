"""
Constructs a *Tanner graph* from a given parity-check matrix `H`, where rows
correspond to *check nodes* and columns to *variable nodes*.

The resulting *bipartite* graph indexes variable nodes as `1:n` and check
nodes as `n+1:n+m` for an `m × n` matrix `H`."""
function tanner_graph_from_parity_matrix(H::AbstractMatrix)
    m, n = size(H)
    g = SimpleGraph(n + m)
    var_nodes = collect(1:n)
    check_nodes = collect(n+1:n+m)
    for row in 1:m, col in 1:n
        H[row, col] && add_edge!(g, col, n + row)
    end
    return (graph=g, left=var_nodes, right=check_nodes)
end

"""Reconstructs the parity-check matrix from a Tanner graph `g`.

!!! note
    The first block of vertices corresponds to variable nodes and the remaining to check nodes. It
    returns a parity check matrix `H` of size `(number of check nodes) × (number of variable nodes).`
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

"""
Constructs a bipartite Tanner graph representing the cycle code of length `n`. The graph
consists of `n` variable nodes and `n` check nodes, totaling `2n` nodes. Each check node `i`
(for `i` in `1...n`) connects to variable nodes `i` and (``i+1 \\mod n``). This structure results
in a parity-check matrix `H`, where each row contains exactly two `1`s, encoding the edges of a
length-`n` cycle.
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
Represents the CSS quantum code `Q(G₁ × G₂)` constructed from two binary codes with
parity-check matrices `H₁` and `H₂`, using the hypergraph product formulation introduced
by [tillich2013quantum](@cite).

This construction corresponds to a specific product of Tanner graphs:
- Let `G₁ = T(V₁, C₁, E₁)` and `G₂ = T(V₂, C₂, E₂)` be the Tanner graphs of `H₁` and `H₂`.
- The product graph `G₁ × G₂` has vertex set `(V₁ × V₂) ∪ (C₁ × C₂)` and check set `(C₁ × V₂) ∪ (V₁ × C₂)`.
- The Tanner subgraphs `G₁ ×ₓ G₂` and `G₁ ×𝓏 G₂` define classical codes `Cₓ` and `C𝓏` used in the CSS construction.

The `hgp(H₁, H₂)` function algebraically realizes this graph-theoretic product using Kronecker operations,
yielding the `X`- and `Z`-type parity-check matrices:

- `H_X = [H₁ ⊗ I  |  I ⊗ H₂ᵗ]` corresponds to `G₁ ×ₓ G₂`
- `H_Z = [I ⊗ H₂  |  H₁ᵗ ⊗ I]` corresponds to `G₁ ×𝓏 G₂`

These matrices ensure `H_X * H_Zᵗ = 0`, satisfying the CSS condition.

See: [tillich2013quantum](@cite), Section 4.3 — “The hypergraph connection, product codes”

# 𝑄(𝐺₁ × 𝐺₂)

The `𝑄(𝐺₁ × 𝐺₂)` quantum LDPC codes represent a broader generalization of **quantum expander codes**
which are derived from the Leverrier-Tillich-Zémor construction [tillich2013quantum](@cite).

```jldoctest examples
julia> using QuantumClifford; using QuantumClifford.ECC; using QECCore

julia> using SparseArrays; # hide

julia> H1 = [1 0 1 0; 0 1 0 1; 1 1 0 0];

julia> H2 = [1 1 0;0 1 1];

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
julia> H = parity_matrix(RepCode(3));

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
struct QuantumTannerGraphProduct <: AbstractCSSCode
    """The classical seed code for the the quantum tanner graph code."""
    H1::AbstractMatrix
    """The classical seed code for the the quantum tanner graph code."""
    H2::AbstractMatrix
end

parity_matrix_xz(c::QuantumTannerGraphProduct) = hgp(c.H1, c.H2)

""" Constructs a 𝑄(𝐺₁ × 𝐺₂) quantum Tanner graph product code using cyclic Tanner graphs of length `2m`.

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC;

julia> m = 10;

julia> c = parity_checks(CyclicQuantumTannerGraphProduct(m));

julia> code_n(c), code_k(c)
(800, 2)
```
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
    H1 = parity_matrix_from_tanner_graph(G1.graph, G1.left, G1.right)
    H2 = parity_matrix_from_tanner_graph(G2.graph, G2.left, G2.right)
    return hgp(H1, H2)
end

parity_matrix_x(c::Union{QuantumTannerGraphProduct,CyclicQuantumTannerGraphProduct}) = parity_matrix_xz(c)[1]

parity_matrix_z(c::Union{QuantumTannerGraphProduct,CyclicQuantumTannerGraphProduct}) = parity_matrix_xz(c)[2]
