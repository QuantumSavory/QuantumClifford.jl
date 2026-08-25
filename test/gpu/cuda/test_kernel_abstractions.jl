@testitem "CUDA" tags = [:cuda] begin

    include("../common/test_platform.jl")

    using CUDA: CuArray, synchronize
    const AT = CuArray

    test_platform(AT, synchronize)

end
