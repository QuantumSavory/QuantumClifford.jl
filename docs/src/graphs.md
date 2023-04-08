# [Graph States](@id Graph-States)

!!! warning "The `graphstate` API is not considered stable"
    `graphstate` returns a lot of information about encoding a given stabilizer state in a graph. A different API is being designed that streamlines the work with graph states.

Conversion to and from [graph states](https://en.wikipedia.org/wiki/Graph_state) is possible.

Consider a GHZ state:

```@example
using QuantumClifford # hide
ghz(4)
```

It can be converted to a graph state with [`graphstate`](@ref)

```julia
graphstate(ghz(4))[1]
```

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, GraphMakie, CairoMakie;
f = Figure(resolution=(200,200))
a = Axis(f[1,1])
graphplot!(a,graphstate(ghz(4))[1])
hidedecorations!(a); hidespines!(a)
a.aspect = DataAspect()
save("ghz4graph.png", f); nothing
```

![](ghz4graph.png)

Notice that the initial GHZ state was not in the typical graph state form. We can see that explicitly by converting back and forth between the two forms:

```jldoctest graph
julia> using Graphs, QuantumClifford

julia> ghz(4)
+ XXXX
+ ZZ__
+ _ZZ_
+ __ZZ

julia> Stabilizer(Graph(ghz(4)))
+ XZZZ
+ ZX__
+ Z_X_
+ Z__X
```

There is a set of single-qubit operations that can convert any stabilizer tableau into a state representable as a graph. These transformations are performed implicitly by the `Graph` constructor when converting from a `Stabilizer`. If you need the explicit transformation you can use the [`graphstate`](@ref) function that specifies which qubits require a Hadamard, Inverse Phase, or Phase Flip gate. The [`graph_gatesequence`](@ref) or [`graph_gate`](@ref) helper functions can be used to generate the exact operations:

```jldoctest graph
julia> s = ghz(4)
+ XXXX
+ ZZ__
+ _ZZ_
+ __ZZ

julia> g, h_idx, ip_idx, z_idx = graphstate(s);

julia> gate = graph_gate(h_idx, ip_idx, z_idx, nqubits(s))
X₁ ⟼ + X___
X₂ ⟼ + _Z__
X₃ ⟼ + __Z_
X₄ ⟼ + ___Z
Z₁ ⟼ + Z___
Z₂ ⟼ + _X__
Z₃ ⟼ + __X_
Z₄ ⟼ + ___X

julia> canonicalize!(apply!(s,gate)) == canonicalize!(Stabilizer(g))
true
```

These converters also provides for a convenient way to create graph and cluster states, by using the helper constructors provided in `Graphs.jl`.

```jldoctest graph
julia> Stabilizer(grid([4,1])) # Linear cluster state
+ XZ__
+ ZXZ_
+ _ZXZ
+ __ZX

julia> Stabilizer(grid([2,2])) # Small 2D cluster state
+ XZZ_
+ ZX_Z
+ Z_XZ
+ _ZZX
```

Graphs are represented with the `Graphs.jl` package and plotting can be done both in `Plots.jl` and `Makie.jl` (with `GraphMakie`).