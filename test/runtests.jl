using QuantumClifford
using TestItemRunner
using Pkg

const GROUP = get(ENV, "GROUP", "QuantumClifford")

if GROUP == "QuantumClifford"
    include("runtests_QuantumClifford.jl")
elseif GROUP == "QECCore"
    Pkg.activate(joinpath(dirname(@__DIR__), "lib", GROUP))
    Pkg.test(GROUP)
end