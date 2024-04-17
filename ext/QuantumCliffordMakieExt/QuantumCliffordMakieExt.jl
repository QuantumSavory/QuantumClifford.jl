module QuantumCliffordMakieExt

using Makie
using QuantumClifford
import QuantumClifford: stabilizerplot, stabilizerplot_axis

# If you want to directly use heatmap
function Makie.convert_arguments(P::Type{<:Makie.Heatmap}, s::Stabilizer)
    h = stab_to_gf2(s)
    r = h[:,1:end÷2] + h[:,end÷2+1:end]*2
    r = r[end:-1:1,:]'
    Makie.convert_arguments(P, r)
end

# A complete Makie recipe
Makie.@recipe(
function (scene)
    Makie.Theme(;
        xzcomponents = :together,
        colormap = Makie.cgrad([:lightgray,Makie.RGB(0x1b9e77),Makie.RGB(0xd95f02),Makie.RGB(0x7570b3)], 4, categorical = true),
        colorrange = (-0.5, 3.5)
    )
end,
StabilizerPlot,
stabilizer
)

function Makie.plot!(myplot::StabilizerPlot)
    s = myplot[:stabilizer][]
    r = if myplot[:xzcomponents][]==:split
        QuantumClifford.stab_to_gf2(s)
    elseif myplot[:xzcomponents][]==:together
        h = QuantumClifford.stab_to_gf2(s)
        h[:,1:end÷2] + h[:,end÷2+1:end]*2
    else
        throw(ErrorException("`xzcomponents` should be `:split` or `:together`"))
    end
    r = r[end:-1:1,:]'
    hm = Makie.heatmap!(myplot, r;
        colorrange = (0, 3),
        colormap=myplot.colormap
    )
    for k in [:colorscale, :highclip, :lowclip]
        myplot[k] = hm[k]
    end
    myplot
end

"""Create a complete Makie figure of the tableaux.

This function is a temporary fix for Makie limitations. It lets you make a complete figure.

See [Makie#379](https://github.com/JuliaPlots/Makie.jl/issues/379)."""
function stabilizerplot_axis(subfig, s; colorbar=true, args...)
    ax = Makie.Axis(subfig[1,1])
    p = stabilizerplot!(ax,s; args...)
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    ax.aspect = Makie.DataAspect()
    colorbar && Makie.Colorbar(subfig[1, 2], p, ticks = ((0.5:3.51)*3/4, ["I", "X", "Z", "Y"]), vertical = true, flipaxis = true)
    #Makie.colsize!(subfig.layout, 1, Makie.Aspect(1, min(1,size(s,2)/size(s,1))))
    subfig,ax,p
end

end
