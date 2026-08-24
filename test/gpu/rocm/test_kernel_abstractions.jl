@testitem "ROCm" tags = [:rocm] begin

    include("../common/test_platform.jl")

    using AMDGPU: ROCArray, devices, synchronize
    const AT = ROCArray

    const can_run = length(devices()) > 0

    @testset "Device availability" begin
        @test can_run
    end

    if can_run
        test_platform(AT, synchronize)
    end

end
