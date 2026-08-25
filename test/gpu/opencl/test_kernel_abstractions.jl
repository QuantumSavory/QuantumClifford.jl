@testitem "OpenCL" tags = [:opencl] begin

    include("../common/test_platform.jl")

    import pocl_jll
    using OpenCL: CLArray, cl.finish, cl.queue
    const AT = CLArray

    # TODO: Revisit this once the POCL code generation issues are resolved.
    block_sizes = fill(256, round_count)
    synchronize() = finish(queue())
    test_platform(AT, synchronize)

end
