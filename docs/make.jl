push!(LOAD_PATH,"../src/")

using Revise # for interactive doc updates
using Documenter
using DocumenterCitations
using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits
using QuantumInterface
using QECCore
using QuantumClifford.ECC

ENV["HECKE_PRINT_BANNER"] = "false"
import Hecke
const QuantumCliffordHeckeExt = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)

ENV["OSCAR_PRINT_BANNER"] = "false"
import Oscar
const QuantumCliffordOscarExt = Base.get_extension(QuantumClifford, :QuantumCliffordOscarExt)

import JuMP
const QuantumCliffordJuMPExt = Base.get_extension(QuantumClifford, :QuantumCliffordJuMPExt)

#DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)

ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
ENV["COLUMNS"] = 80

bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"),style=:authoryear)

makedocs(
plugins = [bib],
doctest = false,
clean = true,
sitename = "QuantumClifford.jl",
format = Documenter.HTML(size_threshold_ignore = ["API.md", "ECC_API.md"]),
modules = [QuantumClifford, QuantumClifford.Experimental.NoisyCircuits, QuantumClifford.ECC, QuantumInterface, QuantumCliffordHeckeExt, QuantumCliffordOscarExt, QuantumCliffordJuMPExt, QECCore],
warnonly = [:missing_docs],
linkcheck = true,
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
"Circuit Simulation" => [
    "Simulation of Noisy Circuits" => "noisycircuits.md",
    "Monte Carlo" => "noisycircuits_mc.md",
    "Perturbative Expansions" => "noisycircuits_perturb.md",
    "ECC example" => "ecc_example_sim.md",
    "Circuit Operations" => "noisycircuits_ops.md",
    "API" => "noisycircuits_API.md"
],
"ECC compendium" => [
    "Evaluating codes and decoders" => "ECC_evaluating.md",
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
