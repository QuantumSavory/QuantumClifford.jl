@testitem "OpenCL" tags = [:opencl] begin

    include("implementation/test_platform.jl")

    import pocl_jll
    using OpenCL: CLArray, cl.devices, cl.platforms, cl.finish, cl.queue
    const AT = CLArray

    const can_run = any(
        length(devices(platform)) > 0 for platform in platforms()
            )

    @testset "Device availability" begin
        @test can_run
    end

    if can_run
        synchronize() = finish(queue())
        test_platform(AT, synchronize)
    end

end
