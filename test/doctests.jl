using Documenter

function doctests()
    if VERSION >= v"1.7"
        @testset "Doctests" begin
            DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)
            doctest(QuantumClifford)
        end
    end
end

doctests()