using Documenter
using QuantumClifford

ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
ENV["COLUMNS"] = 80
@testset "Doctests" begin
    DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford; using QuantumClifford.ECC); recursive=true)
    doctest(QuantumClifford,
    #fix=true
    )
end
