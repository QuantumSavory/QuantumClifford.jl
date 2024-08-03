@testitem "Doctests" tags=[:doctests] begin
    using Documenter
    using QuantumClifford

    ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
    ENV["COLUMNS"] = 80
    DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford; using QuantumClifford.ECC); recursive=true)
    doctestfilters = [r"(QuantumClifford\.|)"]
    doctest(QuantumClifford;
            doctestfilters
            #fix=true
           )
    # TODO failures in VSCode related to https://github.com/julia-vscode/TestItemRunner.jl/issues/49
end
