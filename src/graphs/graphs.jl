import Graphs: Graphs, nv, AbstractGraph

const clifford_id1 = one(CliffordOperator, 1)

# elided as it's a non-public API and might introduce unexpected behavior for downstream users
# we didn't wrap it inside a new struct as that introduces boiler plate for show etc.
function _apply_vop!(vops::Vector{SingleQubitOperator}, op::AbstractSingleQubitOperator; phases::Bool=true)
    q = affectedqubits(op)[1]
    # Convert into a single qubit operator acting on index 1
    op_prime = SingleQubitOperator(SingleQubitOperator(op), 1)
    # A trick that performs slightly better than directly op_prime * CliffordOperator(vops[q])
    # as converting to `CliffordOperator` is expensive
    vops[q] = SingleQubitOperator(op_prime * (vops[q] * clifford_id1))
end

# NOTE: Not tested yet, so commented out for now
# function _apply_vop_right!(vops::Vector{SingleQubitOperator}, op::AbstractSingleQubitOperator; phases::Bool=true)
#     q = affectedqubits(op)[1]
#     # Convert into a single qubit operator acting on index 1
#     op_prime = SingleQubitOperator(SingleQubitOperator(op), 1)
#     vops[q] = SingleQubitOperator(vops[q] * (op_prime * clifford_id1))
# end

# TODO: Implement sanity_check like the one in Stim
struct GraphState{G, S}
    # Graph representing the graph state
    graph::G
    # Vertex Operators that need to be applied to a stabilizer state
    # before the state is representable as a "pure" graph state
    vops::Vector{S}
end

nqubits(g::GraphState) = nv(g.graph)

"""Return the underlying graph of the graph state"""
graph(g::GraphState) = g.graph

"""Return the VOPs (Vertex Operators) of the graph state"""
vops(g::GraphState) = g.vops

function GraphState(graph::G, h_idx::Vector{Int}, ip_idx::Vector{Int}, z_idx::Vector{Int}) where G
    vops = [SingleQubitOperator(sId1(1)) for _ in 1:nv(graph)]
    for id in h_idx
        _apply_vop!(vops, sHadamard(id))
    end
    for id in ip_idx
        _apply_vop!(vops, sInvPhase(id))
    end
    for id in z_idx
        _apply_vop!(vops, sZ(id))
    end
    GraphState{G, SingleQubitOperator}(graph, vops)
end

function GraphState(stab::Stabilizer)
    GraphState(graphstate(stab)...)
end

function Stabilizer(self::GraphState{<:AbstractGraph})
    s_raw = Stabilizer(self.graph)
    for id in eachindex(self.vops)
        apply!(s_raw, SingleQubitOperator(inv(self.vops[id]), id))
    end
    s_raw
end

"""An in-place version of [`graphstate`](@ref)."""
function graphstate!(stab::Stabilizer)
    n = nqubits(stab)
    # canonicalization is in-place
    _, r, _, permx, permz = canonicalize_gott!(stab)
    # This is equivalent to saying perm(i) = (permz ∘ permx)(i)
    # No typo above. permz ∘ permx because in Gottesman canonicalization
    # we first do Gaussian elimination on the X half and then on the Z half
    # when writing tableau in the matrix format.
    perm = permx[permz]

    # Swap X ↔ Z in the (r+1):n columns so we have only X,Y,I diagonal and only Z off-diagonal
    h_idx = [perm[i] for i in (r+1):n]
    # Swap Y → X on the diagonal so we have only X, I on the diagonal
    ip_idx = [perm[i] for i in 1:n if stab[i,i]==(true,true)]
    # Since ZXZ = -X, apply Z on rows with negative phase to make all phases positive.
    phase_flips = [perm[i] for i in 1:n if phases(stab)[i]!=0x0]

    # canonicalized stabilizer tableau serves as an adjacency matrix
    graph = Graphs.SimpleGraphFromIterator((
        Graphs.SimpleEdge(perm[i],perm[j])
        for i in 1:n, j in 1:n
        if i!=j && stab[i,j]!=(false,false)
    ))
    Graphs.add_vertices!(graph, n-Graphs.nv(graph))
    graph, h_idx, ip_idx, phase_flips
end

""" Convert any stabilizer state to a graph state

[Graph states](https://en.wikipedia.org/wiki/Graph_state) are a special type
of entangled stabilizer states that can be represented by a graph.
For a graph ``G=(V,E)`` the corresponding stabilizers are ``S_v = X_v \\prod_{u ∈ N(v)} Z_u``.
Notice that such tableau rows contain only a single X operator.
There is a set of single qubit gates that converts any stabilizer state to a graph state.

This function returns the graph state corresponding to a stabilizer and the gates that
might be necessary to convert the stabilizer into a state representable as a graph.

For a tableau `stab` you can convert it with:
```julia
graph, hadamard_idx, iphase_idx, flips_idx = graphstate()
```
where `graph` is the graph representation of `stab`,
and the rest specifies the single-qubit gates converting `stab` to `graph`:
`hadamard_idx` are the qubits that require a Hadamard gate (mapping X ↔ Z),
`iphase_idx` are (different) qubits that require an inverse Phase gate (Y → X),
and `flips_idx` are the qubits that require a phase flip (Pauli Z gate),
after the previous two sets of gates.

```jldoctest graph
julia> using Graphs

julia> s = S" XXX
              ZZ_
             -_ZZ";


julia> g, h_idx, ip_idx, z_idx = graphstate(s);


julia> collect(edges(g))
2-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
 Edge 1 => 2
 Edge 1 => 3

julia> h_idx
2-element Vector{Int64}:
 2
 3

julia> ip_idx
Int64[]

julia> z_idx
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
performs transformations in-place, use `graphstate!`. It would
perform `canonicalize_gott!` on its argument as it finds a way
to convert it to a graph state.
"""
graphstate(s::AbstractStabilizer) = graphstate!(copy(stabilizerview(s)))

Graphs.Graph(s::AbstractStabilizer) = graphstate(s)[1]

"""Convert a graph representing a stabilizer state to an explicit Stabilizer.

See also: [`graphstate`](@ref)"""
function Stabilizer(g::Graphs.Graph)
    s = zero(Stabilizer, Graphs.nv(g)) # number of vertices
    for v in Graphs.vertices(g)
        s[v,v]=(true,false)
        for n in Graphs.neighbors(g,v)
            s[v,n]=(false,true)
        end
    end
    return s
end

"""A helper function converting the gate indices from [`graphstate`](@ref) into a sequence of gates.

```jldoctest
julia> s = S" XXX
              YZ_
             -_ZZ";


julia> graph, h_idx, ip_idx, z_idx = graphstate(s);


julia> gates = graph_gatesequence(h_idx, ip_idx, z_idx);


julia> for gate in vcat(gates...) apply!(s, gate) end


julia> s # This is now a graph state (notice you need to multiply row 1 by row 2)
+ YYZ
+ XZ_
+ _ZX

julia> canonicalize!(s) == canonicalize!(Stabilizer(graph))
true
```

See also: [`graph_gatesequence`](@ref)
"""
function graph_gatesequence(h_idx::Vector{Int}, ip_idx::Vector{Int}, z_idx::Vector{Int})
    ([sHadamard(i) for i in h_idx], [sInvPhase(i) for i in ip_idx], [sZ(i) for i in z_idx])
end

"""A helper function converting the gate indices from [`graphstate`](@ref) into a Clifford operator.

```jldoctest
julia> s = S" XXX
              YZ_
             -_ZZ";


julia> graph, h_idx, ip_idx, z_idx = graphstate(s);


julia> gate = graph_gate(h_idx, ip_idx, z_idx, nqubits(s));


julia> apply!(s, gate) # This is now a graph state (notice you need to multiply row 1 by row 2)
+ YYZ
+ XZ_
+ _ZX

julia> canonicalize!(s) == canonicalize!(Stabilizer(graph))
true
```

See also: [`graph_gatesequence`](@ref)
"""
function graph_gate(h_idx, ip_idx, z_idx, n)
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
