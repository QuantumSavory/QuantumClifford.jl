import RecipesBase

RecipesBase.@recipe function f(s::Stabilizer; xzcomponents=:split)
    seriestype  := :heatmap
    aspect_ratio := :equal
    yflip := true
    colorbar := false
    grid := false
    framestyle := :none
    if xzcomponents==:split
        stab_to_gf2(s)
    elseif xzcomponents==:together
        h = stab_to_gf2(s)
        h[:,1:end÷2]*2 + h[:,end÷2+1:end]
    else
        throw(ErrorException("`xzcomponents` should be `:split` or `:together`"))
    end
end
