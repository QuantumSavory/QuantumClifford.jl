push!(LOAD_PATH,"../src/")

using Documenter
using SimpleClifford

DocMeta.setdocmeta!(SimpleClifford, :DocTestSetup, :(using SimpleClifford); recursive=true)

makedocs(
    sitename = "SimpleClifford",
    format = Documenter.HTML(),
    modules = [SimpleClifford]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
