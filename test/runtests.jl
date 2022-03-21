using Test, Random
using QuantumClifford
using QuantumClifford: stab_looks_good, mixed_stab_looks_good, destab_looks_good, mixed_destab_looks_good
using QuantumClifford: apply_single_x!, apply_single_y!, apply_single_z!
using QuantumClifford: mul_left!
using QuantumClifford: perm_inverse
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

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

doset("paulis")             && include("./test_paulis.jl")
doset("stabilizers")        && include("./test_stabs.jl")
doset("canonicalization")   && include("./test_stabcanon.jl")
doset("inner")              && include("./test_inner.jl")
doset("gf2")                && include("./test_gf2.jl")
doset("projection")         && include("./test_projections.jl")
doset("expect")             && include("./test_expect.jl")
doset("trace")              && include("./test_trace.jl")
doset("cliffords")          && include("./test_cliff.jl")
doset("symbolic cliffords") && include("./test_symcliff.jl")
doset("random")             && include("./test_random.jl")
doset("noisy circuits")     && include("./test_noisycircuits.jl")
doset("allocations")        && include("./test_allocations.jl")
doset("bitpack")            && include("./test_bitpack.jl")
doset("graphs")             && include("./test_graphs.jl")
doset("doctests")           && include("./doctests.jl")

using Aqua
doset("aqua") && Aqua.test_all(QuantumClifford, ambiguities=false,project_toml_formatting=false)
