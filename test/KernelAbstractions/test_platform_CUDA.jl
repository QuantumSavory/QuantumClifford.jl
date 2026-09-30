
#=============================================================================#
@testitem "CUDA" tags = [:cuda] begin

    include("implementation/test_platform.jl")

    using CUDA: CuArray, devices, synchronize
    const AT = CuArray

    const can_run = length(devices()) > 0

    @testset "Device availability" begin
        @test can_run
    end

    if can_run
        test_platform(synchronize, AT)
    end

end
#=============================================================================#
