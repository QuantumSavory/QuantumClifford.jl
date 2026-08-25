@testitem "ROCm" tags = [:rocm] begin

    include("../common/test_platform.jl")

    using AMDGPU: ROCArray, synchronize
    const AT = ROCArray

    test_platform(AT, synchronize)

end
