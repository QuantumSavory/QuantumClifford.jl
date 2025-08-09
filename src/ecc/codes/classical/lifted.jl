"""Classical codes lifted over a group algebra, used for lifted product code construction ([panteleev2021degenerate](@cite), [panteleev2022asymptotically](@cite))

Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function LiftedCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw(THROW_MISSING_HECKE)
    end
    return ext.LiftedCode(args...; kwargs...)
end
