include("imports.jl")
include("definitions.jl")
include("utilities.jl")
include("benchmark_KA_mul_leftright.jl")
include("benchmark_KA_canonicalize.jl")

@inline function benchmark_platform(AT, synchronize, path)
    benchmark_KA_mul_leftright(AT, synchronize, path; phases = Val(true))
    benchmark_KA_mul_leftright(AT, synchronize, path; phases = Val(false))
    benchmark_KA_canonicalize_all(AT, synchronize, path)
end
