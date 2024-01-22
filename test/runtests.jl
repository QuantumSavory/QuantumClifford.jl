using SafeTestsets
using QuantumClifford
using DotEnv

DotEnv.config(path = ".env")

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
println("ENV[\"PYTHON\"] = \"$(get(ENV,"PYTHON",nothing))\"")

# in order to run the gpu tests automatically set GPU_TESTS to true in the .env file
if ENV["GPU_TESTS"] == "true"
    @doset "gpu"
else
    println("skipping gpu tests (set GPU_TESTS=true to test gpu)")
end

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
@doset "ecc_decoder_all_setups"
@doset "ecc_encoding"
@doset "ecc_syndromes"
@doset "precompile"
@doset "pauliframe"
@doset "allocations"
VERSION >= v"1.10" && @doset "doctests"
get(ENV,"JET_TEST","")=="true" && @doset "jet"
VERSION >= v"1.10" && @doset "aqua"
