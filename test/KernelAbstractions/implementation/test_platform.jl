
#=============================================================================#
include("imports.jl")
include("definitions.jl")
include("utilities.jl")

include("suites/test_KA_mul.jl")
include("suites/test_KA_canonicalization.jl")

@inline function test_platform(synchronize, AT)::Nothing
    cache = AllocCache()

    @testset "mul" begin
        test_KA_mul(synchronize, AT, cache)
    end
    @testset "canonicalization" begin
        test_KA_canonicalization(synchronize, AT, cache)
    end

    return nothing
end
#=============================================================================#
