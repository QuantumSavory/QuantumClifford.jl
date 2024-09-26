"""Classical codes lifted over a group algebra, used for lifted product code construction ([panteleev2021degenerate](@cite), [panteleev2022asymptotically](@cite))

Implemented as a package extension with Hecke. Check the QuantumClifford documentation for more details on that extension."""
function LiftedCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `LiftedCode` depends on the package `Hecke` but you have not installed or imported it yet. Immediately after you import `Hecke`, the `LiftedCode` will be available.")
    end
    return ext.LiftedCode(args...; kwargs...)
end
