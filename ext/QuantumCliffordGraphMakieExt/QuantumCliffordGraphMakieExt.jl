module QuantumCliffordGraphMakieExt

using GraphMakie
using Makie: distinguishable_colors
using GraphMakie: Point2f
using Graphs
using QuantumClifford
using QuantumClifford.ECC

"""Plots the Tanner graphs associated with quantum tanner graph product codes and quantum expander codes.

Requires a Makie.jl backend to be loaded, e.g. `using CairoMakie`.

# Example

Here is an example of quantum tanner graph product code with repetition code as classical seed.

```@example
julia> using GraphMakie # hide

julia> using Makie; using CairoMakie; # hide

julia> using Graphs # hide

julia> using QuantumClifford # hide

julia> using QuantumClifford.ECC # hide

julia> using SparseArrays; # hide

julia> H1 = sparse(parity_checks(RepCode(3)));

julia> H2 = sparse(parity_checks(RepCode(4)));

julia> G1 = tanner_graph_from_parity_matrix(H1);

julia> G2 = tanner_graph_from_parity_matrix(H2);

julia> PG_X = product_tanner_graph_X(G1, G2);

julia> PG_Z = product_tanner_graph_Z(G1, G2);

julia> graph_X = PG_X.graph;

julia> graph_Z = PG_Z.graph;

julia> Graphs.is_bipartite(graph_X)
true

julia> Graphs.is_bipartite(graph_Z)
true

julia> fig, ax, plt = plot_product_tanner_graph(PG_X);

julia> fig

julia> fig1, ax, plt = plot_product_tanner_graph(PG_Z);

julia> fig1
```
"""
function QuantumClifford.ECC.plot_product_tanner_graph(ptg::NamedTuple; filename="my_graph.png")
    g = ptg.graph
    n_vars = length(ptg.var_nodes)
    n_checks = length(ptg.check_nodes)
    total = nv(g)
    var_positions = [Point2f(i, 0.0) for i in 1:n_vars]
    check_positions = [Point2f(i, 1.0) for i in 1:n_checks] 
    positions = vcat(var_positions, check_positions)
    node_shapes = [i ≤ n_vars ? :circle : :rect for i in 1:total]
    var_color = :lightblue
    check_colors = distinguishable_colors(n_checks)
    node_colors = vcat(fill(var_color, n_vars), check_colors) 
    node_stroke_color = :black
    edges = collect(Graphs.edges(g))
    edge_colors = [check_colors[Graphs.dst(e) - n_vars] for e in edges]
    fig, ax, plt = graphplot(g, layout=positions, 
                              node_marker=node_shapes, 
                              node_color=node_colors, 
                              node_strokec=node_stroke_color,
                              node_size=15,
                              edge_color=edge_colors,
                              edge_width=2.0)
    return fig, ax, plt
end

end
