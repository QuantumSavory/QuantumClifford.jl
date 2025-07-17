"""Lifted product codes ([panteleev2021degenerate](@cite), [panteleev2022asymptotically](@cite))

Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function LPCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `LPCode` depends on the package `Hecke` but you have not installed or imported it yet. Immediately after you import `Hecke`, the `LPCode` will be available.")
    end
    return ext.LPCode(args...; kwargs...)
end

"""Implemented in a package extension with Hecke. Check the docs for the [Hecke extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Hecke.jl)"""
function two_block_group_algebra_codes end

"""Implemented in a package extension with Hecke. Check the docs for the [Hecke extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Hecke.jl)"""
function generalized_bicycle_codes end

"""Implemented in a package extension with Hecke. Check the docs for the [Hecke extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Hecke.jl)"""
function bicycle_codes end

"""Implemented in a package extension with Hecke. Check the docs for the [Hecke extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Hecke.jl)"""
function haah_cubic_codes end

"""Implemented in a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function twobga_from_direct_product end

"""Implemented in a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function twobga_from_fp_group end

"""Implemented in a package extension with Hecke."""
function honeycomb_color_codes end

"""Lacross codes ([pecorari2025high](@cite))
Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function Lacross(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `Lacross` quantum LDPC code depends on the package `Hecke` but you have not installed or imported it yet. Immediately after you import `Hecke`, the `Lacross` will be available.")
    end
    return ext.Lacross(args...; kwargs...)
end
