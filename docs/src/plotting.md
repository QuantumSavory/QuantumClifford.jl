# Plotting

Stabilizers have a plot recipe that can be used with `Plots.jl`.

```@repl
using Random; Random.seed!(42); using Plots; using SimpleClifford
plot(random_stabilizer(15,20));
savefig("randstab.png"); nothing
```

```julia
julia> plot(random_stabilizer(15,20));
```

![](randstab.png)
