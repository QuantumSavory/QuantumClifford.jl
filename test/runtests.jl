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

@doset "throws"
@doset "paulis"
@doset "stabs"
@doset "stabcanon"
@doset "mul_leftright"
@doset "inner"
@doset "embed"
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
@doset "memorylayout"
@doset "graphs"
@doset "hash"
@doset "entanglement"
@doset "enumerate"
@doset "quantumoptics"
@doset "ecc"
@doset "precompile"
@doset "pauliframe"
@doset "allocations"
VERSION >= v"1.9" && @doset "doctests"
get(ENV,"JET_TEST","")=="true" && @doset "jet"
VERSION >= v"1.9" && @doset "aqua"
