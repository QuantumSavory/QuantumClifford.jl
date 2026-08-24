@testset "Doctests" begin
    using Documenter
    using QuantumClifford
    using QuantumInterface

    extensions = []

    import Hecke
    push!(extensions, Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt))

    oscar_ext = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)
    if !isnothing(oscar_ext)
        push!(extensions, oscar_ext)
    end

    import JuMP
    push!(extensions, Base.get_extension(QuantumClifford, :QuantumCliffordJuMPExt))

    withenv("LINES" => "80", "COLUMNS" => "80") do
        DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford; using QuantumClifford.ECC); recursive=true)
        modules = [QuantumClifford, QuantumClifford.ECC, QuantumInterface, extensions...]
        doctestfilters = [r"(QuantumClifford\.|)"]
        doctest(nothing, modules;
                doctestfilters
                #fix=true
               )
    end
end
