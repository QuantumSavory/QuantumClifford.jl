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

"""La-cross codes ([pecorari2025high](@cite))
Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function LaCross(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `LaCross` quantum LDPC code depends on the package `Hecke` but you have not installed or imported it yet. Immediately after you import `Hecke`, the `LaCross` will be available.")
    end
    return ext.LaCross(args...; kwargs...)
end

""" Generalized Hypergraph Product codes [panteleev2021degenerate](@cite)
Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function GeneralizedHyperGraphProductCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `GeneralizedHyperGraphProductCode` quantum LDPC code depends on the package `Hecke` but you have not installed or imported it yet. Immediately after you import `Hecke`, the `GeneralizedHyperGraphProductCode` will be available.")
    end
    return ext.GeneralizedHyperGraphProductCode(args...; kwargs...)
end

"""Generalized Bicycle codes ([koukoulekidis2024smallquantumcodesalgebraic](@cite))

Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function GeneralizedBicycleCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `GeneralizedBicycleCode` depends on the package `Hecke` but you have not installed or imported it yet. Immediately after you import `Hecke`, the `GeneralizedBicycleCode` will be available.")
    end
    return ext.GeneralizedBicycleCode(args...; kwargs...)
end

"""Extended Generalized Bicycle codes ([koukoulekidis2024smallquantumcodesalgebraic](@cite))

Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function ExtendedGeneralizedBicycleCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `ExtendedGeneralizedBicycleCode` depends on the package `Hecke` but you have not installed or imported it yet. Immediately after you import `Hecke`, the `ExtendedGeneralizedBicycleCode` will be available.")
    end
    return ext.ExtendedGeneralizedBicycleCode(args...; kwargs...)
end
