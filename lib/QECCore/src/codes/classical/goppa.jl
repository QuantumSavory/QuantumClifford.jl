"""Goppa codes [goppa1970new](@cite)."""
function Goppa(args...; kwargs...)
    ext = Base.get_extension(QECCore, :QECCoreNemoExt)
    if isnothing(ext)
        throw("The `Goppa` depends on the package `Nemo` but you have not installed or imported it yet. Immediately after you import `Nemo`, the `Goppa` will be available.")
    end
    return ext.Goppa(args...; kwargs...)
end

function random_Goppa_code end
