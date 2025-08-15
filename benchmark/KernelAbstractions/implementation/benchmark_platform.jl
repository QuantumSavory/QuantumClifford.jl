include("imports.jl")
include("definitions.jl")
include("utilities.jl")
include("benchmark_KA_mul_leftright.jl")

@inline function benchmark_platform(AT, synchronize, path)
    path *= "/" * string(value(now()) - UNIXEPOCH)
    benchmark_KA_mul_leftright(AT, synchronize, path; phases = Val(true))
    benchmark_KA_mul_leftright(AT, synchronize, path; phases = Val(false))
end
