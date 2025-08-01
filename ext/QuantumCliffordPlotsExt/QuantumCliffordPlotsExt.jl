module QuantumCliffordPlotsExt

using Plots
using QuantumClifford
include("../../src/throws.jl")

@recipe function f(s::QuantumClifford.Stabilizer; xzcomponents=:together)
    seriestype  := :heatmap
    aspect_ratio := :equal
    yflip := true
    colorbar := :none
    colorbar_discrete_values := true
    colorbar_ticks := [0,1,2,3]
    colorbar_formatter := i->["I","X","Z","Y"][Int(floor(i+1))]
    clims := (0,3)
    color_palette := [:blue, :green,:red,:black]
    grid := false
    framestyle := :none
    if xzcomponents==:split
        QuantumClifford.stab_to_gf2(s)
    elseif xzcomponents==:together
        h = QuantumClifford.stab_to_gf2(s)
        h[:,1:end÷2] + h[:,end÷2+1:end]*2
    else
        throw(ErrorException(THROW_INVALID_XZ_COMPONENTS))
    end
end

end
