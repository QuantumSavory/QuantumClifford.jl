"""Lifted product codes ([panteleev2021degenerate](@cite), [panteleev2022asymptotically](@cite))

Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function LPCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw(THROW_MISSING_PACKAGE("LPCode", "Hecke"))
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

"""La-cross codes ([pecorari2025high](@cite))
Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function LaCross(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw(THROW_MISSING_PACKAGE("LaCross", "Hecke"))
    end
    return ext.LaCross(args...; kwargs...)
end
