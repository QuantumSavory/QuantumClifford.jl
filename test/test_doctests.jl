@testitem "Doctests" tags=[:doctests] begin
    using Documenter
    using QuantumClifford
    using QuantumClifford.Experimental.NoisyCircuits
    using QuantumInterface

    ENV["HECKE_PRINT_BANNER"] = "false"
    import Hecke
    const QuantumCliffordHeckeExt = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)

    ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
    ENV["COLUMNS"] = 80
    DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford; using QuantumClifford.ECC); recursive=true)
    modules = [QuantumClifford, QuantumClifford.Experimental.NoisyCircuits, QuantumClifford.ECC, QuantumInterface, QuantumCliffordHeckeExt]
    doctestfilters = [r"(QuantumClifford\.|)"]
    doctest(nothing, modules;
            doctestfilters
            #fix=true
           )
    # TODO failures in VSCode related to https://github.com/julia-vscode/TestItemRunner.jl/issues/49
end
