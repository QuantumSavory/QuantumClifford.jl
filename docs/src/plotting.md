# Plotting

Stabilizers have a plot recipe that can be used with `Plots.jl`. It simply displays the corresponding parity check matrix (extracted with [`stab_to_gf2`](@ref)) as a bitmap image.

```@eval
using Random; Random.seed!(1); using Plots; using QuantumClifford
plot(random_stabilizer(40,50));
savefig("plot-randstab.png"); nothing
```

```julia
julia> plot(random_stabilizer(40,50));
```

![](plot-randstab.png)

```@eval
using Random; Random.seed!(1); using Plots; using QuantumClifford
plot(canonicalize!(random_stabilizer(40,50)));
savefig("plot-canostab.png"); nothing
```

```julia
julia> plot(canonicalize!(random_stabilizer(40,50)));
```

![](plot-canostab.png)

```@eval
using Random; Random.seed!(1); using Plots; using QuantumClifford
plot(canonicalize_gott!(random_stabilizer(40,50))[1]);
savefig("plot-gottstab.png"); nothing
```

```julia
julia> plot(canonicalize_gott!(random_stabilizer(40,50))[1]);
```

![](plot-gottstab.png)


```@eval
using Random; Random.seed!(1); using Plots; using QuantumClifford
plot(canonicalize_gott!(random_stabilizer(40,50))[1]; xzcomponents=:together);
savefig("plot-gottstab-together.png"); nothing
```

```julia
julia> plot(canonicalize_gott!(random_stabilizer(40,50))[1]; xzcomponents=:together);
```

![](plot-gottstab-together.png)

```@eval
using Random; Random.seed!(1); using Plots; using QuantumClifford
plot(canonicalize_rref!(random_stabilizer(40,50),1:50)[1]; xzcomponents=:together);
savefig("plot-rref-together.png"); nothing
```

```julia
julia> plot(canonicalize_rref!(random_stabilizer(40,50),1:50)[1]; xzcomponents=:together);
```

![](plot-rref-together.png)
