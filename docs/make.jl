push!(LOAD_PATH,"../src/")

using Documenter
using QuantumClifford

DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)

makedocs(
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
"Datastructure Choice" => "datastructures.md",
"Useful States and Operators" => "commonstates.md",
"Plotting" => "plotting.md",
"API" => "API.md"
]
)

deploydocs(
    repo = "github.com/Krastanov/QuantumClifford.git"
)
