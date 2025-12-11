
#=============================================================================#
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
        # TODO: Revisit this once the POCL code generation issues are resolved.
        block_sizes = fill(KAExt.default_block_size, max_rounds)
        synchronize() = finish(queue())
        test_platform(synchronize, AT)
    end

end
#=============================================================================#
