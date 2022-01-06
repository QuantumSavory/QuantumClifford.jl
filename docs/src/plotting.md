# Plotting

Stabilizers have a plot recipe that can be used with `Plots.jl` or `Makie.jl`. It simply displays the corresponding parity check matrix (extracted with [`stab_to_gf2`](@ref)) as a bitmap image.

The recipes are implemented in a separate package `QuantumCliffordPlots.jl`.

## `Plots.jl`

In `Plots.jl` we have a simple recipe `plot(s::Stabilizer; xzcomponents=...)`
where `xzcomponents=:split` plots the tableau heatmap in a wide form, X bits on the left, Z bits on the right;
or `xzcomponents=:together` plots them overlapping, with different colors for I, X, Z, and Y.

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, Plots;
plot(random_stabilizer(20,30), xzcomponents=:split);
savefig("plot-randstab.png"); nothing
```

```julia
julia> using QuantumClifford, QuantumCliffordPlots, Plots
julia> plot(random_stabilizer(40,50), xzcomponents=:split);
```

![](plot-randstab.png)

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, Plots
plot(canonicalize!(random_stabilizer(20,30)));
savefig("plot-canostab.png"); nothing
```

```julia
julia> using QuantumClifford, QuantumCliffordPlots, Plots
julia> plot(canonicalize!(random_stabilizer(20,30)), xzcomponents=:split);
```

![](plot-canostab.png)

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, Plots
plot(canonicalize_gott!(random_stabilizer(30))[1], xzcomponents=:split);
savefig("plot-gottstab.png"); nothing
```

```julia
julia> plot(canonicalize_gott!(random_stabilizer(30))[1], xzcomponents=:split);
```

![](plot-gottstab.png)


```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, Plots
plot(canonicalize_gott!(random_stabilizer(30))[1]; xzcomponents=:together);
savefig("plot-gottstab-together.png"); nothing
```

```julia
julia> plot(canonicalize_gott!(random_stabilizer(30))[1]; xzcomponents=:together);
```

![](plot-gottstab-together.png)

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, Plots
plot(canonicalize_rref!(random_stabilizer(20,30),1:30)[1]; xzcomponents=:together);
savefig("plot-rref-together.png"); nothing
```

```julia
julia> plot(canonicalize_rref!(random_stabilizer(20,30),1:30)[1]; xzcomponents=:together);
```

![](plot-rref-together.png)


## `Makie.jl`

Makie's `heatmap` can be directly called on `Stabilizer`.

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, CairoMakie
s = S"IIXZ
      ZZIZ
      YYIZ
      IIIZ
      ZZXZ"
f, ax, p = heatmap(s)
hidedecorations!(ax); hidespines!(ax); # remove ticks and spines
ax.aspect = DataAspect(); # set a one-to-one aspect ratio
save("plot-makie-heatmap.png", f); nothing
```

```julia
using QuantumClifford, QuantumCliffordPlots, CairoMakie
s = S"IIXZ
      ZZIZ
      YYIZ
      IIIZ
      ZZXZ"
f, ax, p = heatmap(s)
hidedecorations!(ax); hidespines!(ax); # remove ticks and spines
ax.aspect = DataAspect(); # set a one-to-one aspect ratio
```

![](plot-makie-heatmap.png)

A full Makie recipe is available as well (supporting `xzcomponents`)

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, CairoMakie
s = S"IIXZ
      ZZIZ
      YYIZ
      IIIZ
      ZZXZ"
f, ax, p = stabilizerplot(s, xzcomponents=:together)
hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()
save("plot-makie-recipe.png", f); nothing
```

```julia
f, ax, p = stabilizerplot(s, xzcomponents=:together)
hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()
```

![](plot-makie-recipe.png)

You might have noticed, Makie recipes do not let you edit the axes or figure,
rather they only permit you to set the plot content.
Which is why we use `hidedecorations!`, `hidesplines!`, and `DataAspect`
to further modify the plot.

You can easily add colorbars (and change the colormap) as well:

```@eval
using Random; Random.seed!(1); using QuantumClifford, QuantumCliffordPlots, CairoMakie
s = S"IIXZ
      ZZIZ
      YYIZ
      IIIZ
      ZZXZ"
fig = Figure()
ax, p = stabilizerplot(fig[1, 1], s, colormap=cgrad(:heat, 4, categorical = true))
hidedecorations!(ax)
hidespines!(ax)
xlims!(ax, 0.5, size(s,2)+0.5) # otherwise there is padding
ylims!(ax, 0.5, size(s,1)+0.5) # otherwise there is padding
# set the aspect ratio of the plot
ax.aspect = DataAspect()
# set the aspect ratio of the layout
colsize!(fig.layout, 1, Aspect(1, size(s,2)/size(s,1))) 
Colorbar(fig[1, 2], p, ticks = (0:3, ["I", "X", "Z", "Y"]))
save("plot-makie-colorbar1.png", fig); nothing
fig = Figure()
ax, p = stabilizerplot(fig[1, 1], s, colormap=cgrad([:lightgray,RGBf(1,0.4,0.4),RGBf(0.3,1,0.5),RGBf(0.4,0.4,1)], 4, categorical = true))
hidedecorations!(ax)
hidespines!(ax)
xlims!(ax, 0.5, size(s,2)+0.5)
ylims!(ax, 0.5, size(s,1)+0.5)
ax.aspect = DataAspect()
colsize!(fig.layout, 1, Aspect(1, size(s,2)/size(s,1))) 
Colorbar(fig[2, 1], p, ticks = (0:3, ["I", "X", "Z", "Y"]), vertical = false, flipaxis = false)
save("plot-makie-colorbar2.png", fig); nothing
```

```julia
fig = Figure()
ax, p = stabilizerplot(fig[1, 1], s, colormap=cgrad(:heat, 4, categorical = true))
hidedecorations!(ax)
hidespines!(ax)
xlims!(ax, 0.5, size(s,2)+0.5) # otherwise there is padding
ylims!(ax, 0.5, size(s,1)+0.5) # otherwise there is padding
# set the aspect ratio of the plot
ax.aspect = DataAspect()
# set the aspect ratio of the layout
colsize!(fig.layout, 1, Aspect(1, size(s,2)/size(s,1))) 
Colorbar(fig[1, 2], p, ticks = (0:3, ["I", "X", "Z", "Y"]))
```

![](plot-makie-colorbar1.png)

Or set a completely custom set of colors:

```julia
fig = Figure()
ax, p = stabilizerplot(fig[1, 1], s, colormap=cgrad([:lightgray,RGBf(1,0.4,0.4),RGBf(0.3,1,0.5),RGBf(0.4,0.4,1)], 4, categorical = true))
hidedecorations!(ax)
hidespines!(ax)
xlims!(ax, 0.5, size(s,2)+0.5)
ylims!(ax, 0.5, size(s,1)+0.5)
ax.aspect = DataAspect()
colsize!(fig.layout, 1, Aspect(1, size(s,2)/size(s,1))) 
Colorbar(fig[2, 1], p, ticks = (0:3, ["I", "X", "Z", "Y"]), vertical = false, flipaxis = false)
```

![](plot-makie-colorbar2.png)