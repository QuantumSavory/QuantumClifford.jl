"""Goppa codes [goppa1970new](@cite)."""
function GoppaCode(args...; kwargs...)
    ext = Base.get_extension(QECCore, :QECCoreNemoExt)
    if isnothing(ext)
        throw("The `GoppaCode` depends on the package `Nemo` but you have not installed or imported it yet. Immediately after you import `Nemo`, the `GoppaCode` will be available.")
    end
    return ext.GoppaCode(args...; kwargs...)
end
