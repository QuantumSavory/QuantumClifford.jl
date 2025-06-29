import Graphs: Graphs, AbstractGraph, nv

"""
[Graph states](https://en.wikipedia.org/wiki/Graph_state) are a special type
of entangled stabilizer states that can be represented by a graph.
For a graph ``G=(V,E)`` the corresponding stabilizers are ``S_v = X_v \\prod_{u ∈ neighbors(v)} Z_u``.
Notice that such tableau are symmetric and contain only X on diagonal and Z off-diagonal.
There is a set of single qubit gates that converts any stabilizer state to a graph state.
"""
# TODO: implement a graph backend interface so this constructor is generic over graph backend
mutable struct GraphState
    graph::AbstractGraph # underlying graph
    # To transform from stabilizer state to graph state, we apply VOP in the orders of H, InvPhase, Z
    # NOTE: Hadamard gate sequence may be swapped with InvPhase sequences. But Z gate sequence shall not be swapped with the previous two.
    h_idx::Vector{Int} # indices for Hadamard gate to transform a given stabilizer state into graph state
    ip_idx::Vector{Int} # indices for InvPhase gate
    z_idx::Vector{Int} # indices for Z gate
end

nqubits(g::GraphState) = nv(g)

"""An in-place version of [`GraphState`](@ref)."""
function GraphState!(stab::Stabilizer)
    n = nqubits(stab)
    # canonicalization is in-place
    _, r, _, permx, permz = canonicalize_gott!(stab)
    perm = permx[permz]

    # Swap X ↔ Z in the (r+1):n columns so we have only X,Y,I diagonal and only Z off-diagonal
    h_idx = [perm[i] for i in (r+1):n]
    # Swap Y → X on the diagonal so we have only X, I on the diagonal
    ip_idx = [perm[i] for i in 1:n if stab[i,i]==(true,true)]
    # Since ZXZ = -X, apply Z on rows with negative phase to make all phases positive.
    z_idx = [perm[i] for i in 1:n if phases(stab)[i]!=0x0]

    # canonicalized stabilizer tableau serves as an adjacency matrix
    graph = Graphs.SimpleGraphFromIterator((
        Graphs.SimpleEdge(perm[i],perm[j])
        for i in 1:n, j in 1:n
        if i!=j && stab[i,j]!=(false,false)
    ))
    Graphs.add_vertices!(graph, n-Graphs.nv(graph))
    return GraphState(graph, h_idx, ip_idx, z_idx)
end

"""A helper function converting the gate indices from [`GraphState`](@ref) into a sequence of gates.

```jldoctest
julia> s = S" XXX
              YZ_
             -_ZZ";


julia> g = GraphState(s);


julia> for gate in graph_gate_sequence(g) apply!(s, gate) end


julia> s # This is now a graph state (notice you need to multiply row 1 by row 2)
+ YYZ
+ XZ_
+ _ZX

julia> canonicalize!(s) == canonicalize!(Stabilizer(g.graph))
true
```

See also: [`graph_gate_sequence`](@ref)
"""
function graph_gate_sequence(self::GraphState)
    vcat([sHadamard(i)::AbstractSingleQubitOperator for i in self.h_idx],
         [sInvPhase(i)::AbstractSingleQubitOperator for i in self.ip_idx],
         [sZ(i)::AbstractSingleQubitOperator for i in self.z_idx])
end

# Return the stabilizer state corresponding to the graph state
function Stabilizer(self::GraphState)
    s_raw = Stabilizer(self.graph)
    for gate in reverse(graph_gate_sequence(self)) apply!(s_raw, inv(gate)) end
    return s_raw
end

""" Convert any stabilizer state to a graph state

This function returns the graph state corresponding to a stabilizer.

```jldoctest graph
julia> using Graphs

julia> s = S" XXX
              ZZ_
             -_ZZ";


julia> g = GraphState(s);


julia> collect(edges(g.graph))
2-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
 Edge 1 => 2
 Edge 1 => 3

julia> g.h_idx
2-element Vector{Int64}:
 2
 3

julia> g.ip_idx
Int64[]

julia> g.z_idx
1-element Vector{Int64}:
 3
```

The `Graphs.jl` library provides many graph-theory tools and the
`MakieGraphs.jl` library provides plotting utilities for graphs.

You can directly call the graph constructor on a stabilizer,
if you just want the graph and do not care about the Clifford
operation necessary to convert an arbitrary state to a state
representable as a graph:

```jldoctest graph
julia> collect(edges( Graph(bell()) ))
1-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
 Edge 1 => 2
```

For a version that does not copy the stabilizer, but rather
performs transformations in-place, use `GraphState!`. It would
perform `canonicalize_gott!` on its argument as it finds a way
to convert it to a graph state.
"""
GraphState(stab::AbstractStabilizer) = GraphState!(copy(stabilizerview(stab)))

function Graphs.Graph(s::AbstractStabilizer)
    stab = copy(stabilizerview(s))
    n = nqubits(stab)
    _, _, _, permx, permz = canonicalize_gott!(stab)
    perm = permx[permz]

    graph = Graphs.SimpleGraphFromIterator((
        Graphs.SimpleEdge(perm[i],perm[j])
        for i in 1:n, j in 1:n
        if i!=j && stab[i,j]!=(false,false)
    ))
    Graphs.add_vertices!(graph, n-Graphs.nv(graph))
    return graph
end

"""Convert a graph representing a stabilizer state to an explicit Stabilizer.

See also: [`GraphState`](@ref)"""
function Stabilizer(g::AbstractGraph)
    s = zero(Stabilizer, Graphs.nv(g)) # number of vertices
    for v in Graphs.vertices(g)
        s[v,v]=(true,false)
        for n in Graphs.neighbors(g,v)
            s[v,n]=(false,true)
        end
    end
    return s
end


"""A helper function converting the gate sequences from [`GraphState`](@ref) into a Clifford operator.

```jldoctest
julia> s = S" XXX
              YZ_
             -_ZZ";


julia> g = GraphState(s);


julia> gate = graph_gate(g);


julia> apply!(s, gate) # This is now a graph state (notice you need to multiply row 1 by row 2)
+ YYZ
+ XZ_
+ _ZX

julia> canonicalize!(s) == canonicalize!(Stabilizer(g.graph))
true
```

See also: [`graph_gate_sequence`](@ref)
"""
function graph_gate(g::GraphState)
    h_idx, ip_idx, z_idx = g.h_idx, g.ip_idx, g.z_idx
    n = nv(g.graph)

    c = one(CliffordOperator,n)
    t = tab(c)
    x = @view t[1:n]
    z = @view t[n+1:end]
    for i in h_idx
        x[i,i] = (false,true)
        z[i,i] = (true,false)
    end
    for i in ip_idx
        x[i,i] = (true,true)
        x.phases[i] = 0x2
    end
    for i in z_idx
        if x[i,i][1]
            x.phases[i] = (x.phases[i]+0x2)&0x3
        end
        if z[i,i][1]
            z.phases[i] = (z.phases[i]+0x2)&0x3
        end
    end
    c
end
