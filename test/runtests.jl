using SafeTestsets
using QuantumClifford

function doset(descr)
    if length(ARGS) == 0
        return true
    end
    for a in ARGS
        if occursin(lowercase(a), lowercase(descr))
            return true
        end
    end
    return false
end

macro doset(descr)
    quote
        if doset($descr)
            @safetestset $descr begin
                include("test_"*$descr*".jl")
            end
        end
    end
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@doset "paulis"
@doset "stabs"
@doset "stabcanon"
@doset "inner"
@doset "gf2"
@doset "projections"
@doset "expect"
@doset "trace"
@doset "cliff"
@doset "symcliff"
@doset "symcontrolled"
@doset "classicalreg"
@doset "random"
@doset "noisycircuits"
@doset "syndromemeas"
@doset "bitpack"
@doset "graphs"
@doset "hash"
@doset "entanglement"
@doset "enumerate"
@doset "ecc"
@doset "precompile"
@doset "pauliframe"
@doset "allocations"
VERSION >= v"1.9" && @doset "doctests"
get(ENV,"JET_TEST","")=="true" && @doset "jet"

using Aqua
VERSION >= v"1.9" && doset("aqua") && begin
    Aqua.test_all(QuantumClifford, project_toml_formatting=false) # TODO: remove project_toml_formatting=false when Aqua is updated
end
