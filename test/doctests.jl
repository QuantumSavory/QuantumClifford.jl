using Documenter

function doctests()
    @testset "Doctests" begin
        DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)
        doctest(QuantumClifford)
    end
end

doctests()
