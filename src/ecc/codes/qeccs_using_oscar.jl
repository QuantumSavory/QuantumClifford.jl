"""Implemented in a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function twobga_from_direct_product end

"""Implemented in a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function twobga_from_fp_group end

"""Implemented in a package extension with `Oscar`. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function boundary_maps end

"""Implemented in a package extension with `Oscar`. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"""
function max_xy_exponents end

"""D-Dimensional Surface codes ([Berthusen_2024](@cite), [Zeng_2019](@cite))

Implemented as a package extension with `Oscar`. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function DDimensionalSurfaceCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `DDimensionalSurfaceCode` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `DDimensionalSurfaceCode` will be available.")
    end
    return ext.DDimensionalSurfaceCode(args...; kwargs...)
end

"""D-Dimensional Toric codes ([Berthusen_2024](@cite), [Zeng_2019](@cite))

Implemented as a package extension with `Oscar`. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function DDimensionalToricCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `DDimensionalToricCode` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `DDimensionalToricCode` will be available.")
    end
    return ext.DDimensionalToricCode(args...; kwargs...)
end

"""Homological Product codes ([Quintavalle_2021](@cite), [xu2024fastparallelizablelogicalcomputation](@cite)).

Implemented as a package extension with `Oscar`. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function HomologicalProductCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `HomologicalProductCode` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `HomologicalProductCode` will be available.")
    end
    return ext.HomologicalProductCode(args...; kwargs...)
end

"""Double Homological Product codes ([Campbell_2019](@cite)).

Implemented as a package extension with `Oscar`. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function DoubleHomologicalProductCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `DoubleHomologicalProductCode` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `DoubleHomologicalProductCode` will be available.")
    end
    return ext.DoubleHomologicalProductCode(args...; kwargs...)
end

"""Generalized Toric codes ([liang2025generalizedtoriccodestwisted](@cite))

Implemented as a package extension with `Oscar`. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function GeneralizedToricCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `GeneralizedToricCode` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `GeneralizedToricCode` will be available.")
    end
    return ext.GeneralizedToricCode(args...; kwargs...)
end

"""Trivariate Tricycle codes ([jacob2025singleshotdecodingfaulttolerantgates](@cite))
Implemented as a package extension with `Oscar`. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function TrivariateTricycleCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `TrivariateTricycleCode` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `TrivariateTricycleCode` will be available.")
    end
    return ext.TrivariateTricycleCode(args...; kwargs...)
end

"""Bivariate Bicycle codes ([jacob2025singleshotdecodingfaulttolerantgates](@cite))
Implemented as a package extension with `Oscar`. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function BivariateBicycleCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `BivariateBicycleCode` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `BivariateBicycleCode` will be available.")
    end
    return ext.BivariateBicycleCode(args...; kwargs...)
end

"""Multivariate Multicycle codes.
Implemented as a package extension with `Oscar`. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function MultivariateMulticycleCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `MultivariateMulticycleCode` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `MultivariateMulticycleCode` will be available.")
    end
    return ext.MultivariateMulticycleCode(args...; kwargs...)
end
