"""
The method `twobga_from_fp_group` provides functionality of forming group algebra via group presentation, extending the capabilities of 2BGA codes.

Implemented as a package extension with Oscar. Check the [QuantumClifford documentation](http://qc.quantumsavory.org/stable/ECC_API/) for more details on that extension."""
function twobga_from_fp_group(args...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if isnothing(ext)
        throw("The `twobga_from_fp_group` depends on the package `Oscar` but you have not installed or imported it yet. Immediately after you import `Oscar`, the `twobga_from_fp_group` will be available.")
    end
    return ext.twobga_from_fp_group(args...)
end
