"""Implemented in a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function d_dimensional_surface_codes end

"""Implemented in a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function d_dimensional_toric_codes end

"""Implemented in a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function pcms end

"""D-Dimensional codes ([Berthusen_2024](@cite))

Implemented as a package extension with Oscar. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function DDimensionalCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `DDimensional` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `DDimensional` will be available.")
    end
    return ext.DDimensionalCode(args...; kwargs...)
end
