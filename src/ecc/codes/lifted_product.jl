"""Lifted product codes ([panteleev2021degenerate](@cite), [panteleev2022asymptotically](@cite))

Implemented as a package extension with Hecke. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function LPCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `LPCode` depends on the package `Hecke` but you have not installed or imported it yet. Immediately after you import `Hecke`, the `LPCode` will be available.")
    end
    return ext.LPCode(args...; kwargs...)
end

"""Implemented in a package extension with Hecke."""
function two_block_group_algebra_codes end

"""Implemented in a package extension with Hecke."""
function generalized_bicycle_codes end

"""Implemented in a package extension with Hecke."""
function bicycle_codes end

"""Implemented in a package extension with Hecke."""
function check_repr_commutation_relation end
