using Test, Random
using QuantumClifford
using QuantumClifford: stab_looks_good, mixed_stab_looks_good, destab_looks_good, mixed_destab_looks_good
using QuantumClifford: apply_single_x!, apply_single_y!, apply_single_z!
using QuantumClifford: mul_left!
using LinearAlgebra: inv
#using Nemo

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

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
            @testset $descr begin include("test_"*$descr*".jl") end
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
@doset "random"
@doset "noisycircuits"
@doset "syndromemeas"
@doset "bitpack"
@doset "graphs"
@doset "hash"
@doset "entanglement"
@doset "enumerate"
VERSION >= v"1.7" && @doset "allocations"
VERSION == v"1.8" && @doset "doctests"
get(ENV,"QUANTUMCLIFFORD_JET_TEST","")=="true" && @doset "jet"

using Aqua
doset("aqua") && begin
    Aqua.test_all(QuantumClifford, ambiguities=false)
    Aqua.test_ambiguities([QuantumClifford,Core]) # otherwise Base causes false positives
end
