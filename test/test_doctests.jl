@testitem "Doctests" tags=[:doctests] begin
    using Documenter
    using QuantumClifford
    using QuantumInterface

    extensions = []

    ENV["HECKE_PRINT_BANNER"] = "false"
    import Hecke
    const QuantumCliffordHeckeExt = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    push!(extensions, QuantumCliffordHeckeExt)

    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        ENV["OSCAR_PRINT_BANNER"] = "false"
        import Oscar
        const QuantumCliffordOscarExt = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
        push!(extensions, QuantumCliffordOscarExt)
    end

    import JuMP
    const QuantumCliffordJuMPExt = Base.get_extension(QuantumClifford, :QuantumCliffordJuMPExt)
    push!(extensions, QuantumCliffordJuMPExt)

    ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
    ENV["COLUMNS"] = 80
    DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford; using QuantumClifford.ECC); recursive=true)
    modules = [QuantumClifford, QuantumClifford.Experimental.NoisyCircuits, QuantumClifford.ECC, QuantumInterface, extensions...]
    doctestfilters = [r"(QuantumClifford\.|)"]
    doctest(nothing, modules;
            doctestfilters
            #fix=true
           )
    # TODO failures in VSCode related to https://github.com/julia-vscode/TestItemRunner.jl/issues/49
end
