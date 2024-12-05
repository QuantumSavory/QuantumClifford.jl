"""
The method `twobga_from_direct_product` constructs the direct product `G₁ × G₂ × … × Gₙ` of multiple groups, where each `Gᵢ` represents a group in the product.

Implemented as a package extension with Oscar. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function twobga_from_direct_product(args...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `twobga_from_direct_product` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `twobga_from_direct_product` will be available.")
    end
    return ext.twobga_from_direct_product(args...)
end
