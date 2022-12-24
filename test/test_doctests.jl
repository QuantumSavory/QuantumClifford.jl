using Documenter
using QuantumClifford

@testset "Doctests" begin
    DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)
    doctest(QuantumClifford)
end
