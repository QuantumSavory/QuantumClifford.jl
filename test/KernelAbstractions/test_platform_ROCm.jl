
#=============================================================================#
@testitem "ROCm" tags = [:rocm] begin

    include("implementation/test_platform.jl")

    using AMDGPU: ROCArray, devices, synchronize
    const AT = ROCArray

    const can_run = length(devices()) > 0

    @testset "Device availability" begin
        @test can_run
    end

    if can_run
        test_platform(synchronize, AT)
    end

end
#=============================================================================#
