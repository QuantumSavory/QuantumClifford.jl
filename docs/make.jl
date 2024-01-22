push!(LOAD_PATH,"../src/")

using Revise # for interactive doc updates
using Documenter
using DocumenterCitations
using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits
using QuantumInterface

#DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)

ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
ENV["COLUMNS"] = 80

bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"),style=:authoryear)

makedocs(
plugins = [bib],
doctest = false,
clean = true,
sitename = "QuantumClifford.jl",
format = Documenter.HTML(size_threshold_ignore = ["API.md"]),
modules = [QuantumClifford, QuantumClifford.Experimental.NoisyCircuits, QuantumClifford.ECC, QuantumInterface],
warnonly = [:missing_docs],
authors = "Stefan Krastanov",
pages = [
"QuantumClifford.jl" => "index.md",
"Stabilizer Tableau Algebra" => [
    "Manual" => "stab-algebra-manual.md",
    "Canonicalization" => "canonicalization.md",
    "Mixed States" => "mixed.md",
    "Graph States" => "graphs.md",
    "Datastructure Choice" => "datastructures.md",
    "Useful States" => "commonstates.md",
],
"Noisy Circuits" => [
    "Simulation of Noisy Circuits" => "noisycircuits.md",
    "Monte Carlo" => "noisycircuits_mc.md",
    "Perturbative Expansions" => "noisycircuits_perturb.md",
    "ECC example" => "ecc_example_sim.md",
    "Circuit Operations" => "noisycircuits_ops.md",
    "API" => "noisycircuits_API.md"
],
"ECC compendium" => [
    "Evaluating codes and decoders" => "ECC_evaluating.md"
    "API" => "ECC_API.md"
],
"All Gates" => "allops.md",
"Visualizations" => "plotting.md",
"API" => "API.md",
"Tutorials and Publications" => "tutandpub.md",
"Suggested Readings & References" => "references.md",
],
)

deploydocs(
    repo = "github.com/QuantumSavory/QuantumClifford.jl.git"
)
