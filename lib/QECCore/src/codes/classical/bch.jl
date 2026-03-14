"""Goppa codes [goppa1970new](@cite)."""
function BCH(args...; kwargs...)
    ext = Base.get_extension(QECCore, :QECCoreNemoExt)
    if isnothing(ext)
        throw("The `BCH` depends on the package `Nemo` but you have not installed or imported it yet. Immediately after you import `Nemo`, the `BCH` will be available.")
    end
    return ext.BCH(args...; kwargs...)
end
