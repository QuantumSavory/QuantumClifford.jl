include("../implementation/test_platform.jl")

using CUDA: CuArray, devices, synchronize
const AT = CuArray
const can_run = length(devices()) > 0

@testset "CUDA" begin

    @testset "Device availability" begin
        @test can_run
    end

    if can_run
        test_platform(AT, synchronize)
    end

end
