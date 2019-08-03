push!(LOAD_PATH,"../src/")

using Documenter
using SimpleClifford

DocMeta.setdocmeta!(SimpleClifford, :DocTestSetup, :(using SimpleClifford); recursive=true)

makedocs(
doctest = false,
clean = true,
sitename = "SimpleClifford.jl",
format = Documenter.HTML(),
modules = [SimpleClifford],
authors = "Stefan Krastanov",
)

deploydocs(
    repo = "https://github.com/Krastanov/SimpleClifford.git"
)
