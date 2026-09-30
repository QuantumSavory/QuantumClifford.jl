
#=============================================================================#
include("imports.jl")
include("definitions.jl")
include("utilities.jl")

include("suites/benchmark_KA_mul.jl")
include("suites/benchmark_KA_canonicalization.jl")

@inline function benchmark_platform(synchronize, AT, path)::Nothing
    # BenchmarkTools evaluates the setup block at the global scope.
    global cache = AllocCache()
    path *= "/" * string(value(now()) - UNIXEPOCH)

    benchmark_KA_mul(synchronize, AT, path)
    benchmark_KA_canonicalization(synchronize, AT, path)

    return nothing
end
#=============================================================================#
