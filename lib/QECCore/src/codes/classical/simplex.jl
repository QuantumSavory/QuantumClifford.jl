"""Simplex code, the dual of the Hamming code [malcolm2025computing](@cite)."""
function Simplex(args...; kwargs...)
    ext = Base.get_extension(QECCore, :QECCoreNemoExt)
    if isnothing(ext)
        throw("The `Simplex` code depends on the package `Nemo` but you have not installed or imported it yet. Immediately after you import `Nemo`, the `Simplex` will be available.")
    end
    return ext.Simplex(args...; kwargs...)
end
