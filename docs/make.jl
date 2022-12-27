push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations
using QuantumClifford

DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"))

makedocs(
bib,
doctest = false,
clean = true,
sitename = "QuantumClifford.jl",
format = Documenter.HTML(),
modules = [QuantumClifford],
authors = "Stefan Krastanov",
pages = [
"QuantumClifford.jl" => "index.md",
"Manual" => "manual.md",
"Canonicalization" => "canonicalization.md",
"Mixed States" => "mixed.md",
"Graph States" => "graphs.md",
"Datastructure Choice" => "datastructures.md",
"Useful States" => "commonstates.md",
"All Gates" => "allops.md",
"Plotting" => "plotting.md",
"API" => "API.md",
"Tutorials and Publications" => "tutandpub.md",
"Noisy Circuits" => [
    "Simulation of Noisy Circuits" => "noisycircuits.md",
    "Monte Carlo" => "noisycircuits_mc.md",
    "Perturbative Expansions" => "noisycircuits_perturb.md",
    "Circuit Operations" => "noisycircuits_ops.md",
    "API" => "noisycircuits_API.md"
],
"Suggested Readings & References" => "references.md",
]
)

deploydocs(
    repo = "github.com/Krastanov/QuantumClifford.jl.git"
)
