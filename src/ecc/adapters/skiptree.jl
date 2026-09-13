"""
SkipTree algorithm (paper Algorithm 2, Appendix E [swaroop2026universal](@cite)).
Returns `(T, P, :H_R)` with `T · G · P ≡ H_R(n) (mod 2)`, where `G` is the
`m × n` edge-vertex incidence of `graph`. `T` is `(n-1) × m` and `(3, 2)`-sparse;
`P` is the `n × n` vertex permutation.

Builds a spanning tree (Kruskal MST) and runs the `LabelFirst`/`LabelLast`
recursion so consecutive labels are within tree-distance 3 (Theorem 7);
row `l` of `T` is the edges of the tree-path between `label[l]` and `label[l+1]`.

The choice of MST and child iteration order affects which valid `T` matrix
comes out, so this does not match the paper's Zenodo `T_X` byte-for-byte.
"""
function skiptree(graph::SimpleGraph{Int}; root::Int = 1)::SkipTreeOutput
    n = nv(graph)
    m = ne(graph)
    n ≥ 2 || throw(ArgumentError("graph must have ≥ 2 vertices (got $n)"))
    1 ≤ root ≤ n || throw(ArgumentError("root must be in 1:$n (got $root)"))
    is_connected(graph) || throw(ArgumentError("graph must be connected"))

    # Edge index in iteration order of edges(graph).
    edge_index = Dict{Tuple{Int,Int},Int}()
    let i = 1
        for e in edges(graph)
            u, v = src(e), dst(e)
            edge_index[(min(u, v), max(u, v))] = i
            i += 1
        end
    end

    mst_edges = kruskal_mst(graph)
    tree = SimpleGraph(n)
    for e in mst_edges
        add_edge!(tree, src(e), dst(e))
    end
    @assert ne(tree) == n - 1

    label = zeros(Int, n)
    visited = falses(n)
    index = Ref(1)

    function label_first(v::Int, skip::Bool)
        visited[v] = true
        label[index[]] = v
        index[] += 1
        children = sort([nbr for nbr in neighbors(tree, v) if !visited[nbr]])
        n_children = length(children)
        for (k, child) in enumerate(children)
            last_in_gen = (k == n_children)
            if last_in_gen && !skip
                label_first(child, false)
            else
                label_last(child)
            end
        end
        return
    end

    function label_last(v::Int)
        visited[v] = true
        for child in sort(collect(neighbors(tree, v)))
            if !visited[child]
                label_first(child, true)
            end
        end
        label[index[]] = v
        index[] += 1
        return
    end

    label_first(root, false)
    @assert index[] == n + 1
    @assert sort(label) == 1:n

    # P[v, l] = 1 iff label[l] == v.
    I_p = Vector{Int}(undef, n)
    J_p = Vector{Int}(undef, n)
    @inbounds for l in 1:n
        I_p[l] = label[l]
        J_p[l] = l
    end
    P = sparse(I_p, J_p, trues(n), n, n)

    # BFS once to get parent/depth for O(depth(u)+depth(v)) LCA path queries.
    parent = zeros(Int, n)
    depth  = zeros(Int, n)
    parent[root] = 0
    bfs_queue = Int[root]
    visited_bfs = falses(n)
    visited_bfs[root] = true
    while !isempty(bfs_queue)
        v = popfirst!(bfs_queue)
        for nbr in neighbors(tree, v)
            if !visited_bfs[nbr]
                visited_bfs[nbr] = true
                parent[nbr] = v
                depth[nbr] = depth[v] + 1
                push!(bfs_queue, nbr)
            end
        end
    end

    function tree_path_edges(u::Int, v::Int)
        # Walk u, v up to their LCA, collecting edges as sorted pairs.
        es = Tuple{Int,Int}[]
        x, y = u, v
        while depth[x] > depth[y]
            push!(es, (min(x, parent[x]), max(x, parent[x])))
            x = parent[x]
        end
        while depth[y] > depth[x]
            push!(es, (min(y, parent[y]), max(y, parent[y])))
            y = parent[y]
        end
        while x != y
            push!(es, (min(x, parent[x]), max(x, parent[x])))
            push!(es, (min(y, parent[y]), max(y, parent[y])))
            x = parent[x]
            y = parent[y]
        end
        es
    end

    I_t = Int[]; J_t = Int[]
    @inbounds for l in 1:(n - 1)
        u, v = label[l], label[l + 1]
        for e in tree_path_edges(u, v)
            push!(I_t, l)
            push!(J_t, edge_index[e])
        end
    end
    T = sparse(I_t, J_t, trues(length(I_t)), n - 1, m)

    SkipTreeOutput(T, P, :H_R)
end
