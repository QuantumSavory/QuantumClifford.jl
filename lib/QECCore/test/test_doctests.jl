@testitem "Doctests" tags=[:doctests] begin
    using Documenter
    using QECCore
    using QuantumClifford

    extensions = []

    ENV["NEMO_PRINT_BANNER"] = "false"
    import Nemo
    const QECCoreNemoExt = Base.get_extension(QECCore, :QECCoreNemoExt)
    push!(extensions, QECCoreNemoExt)

    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        ENV["OSCAR_PRINT_BANNER"] = "false"
        import Oscar
        const QECCoreOscarExt = Base.get_extension(QECCore, :QECCoreOscarExt)
        push!(extensions, QECCoreOscarExt)
    end

    ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
    ENV["COLUMNS"] = 80
    DocMeta.setdocmeta!(QECCore, :DocTestSetup, :(using QECCore; using QuantumClifford); recursive=true)
    modules = [QECCore, extensions...]
    doctestfilters = [r"(QECCore\.|)"]
    doctest(nothing, modules;
            doctestfilters
            #fix=true
           )
    # TODO failures in VSCode related to https://github.com/julia-vscode/TestItemRunner.jl/issues/49
end
