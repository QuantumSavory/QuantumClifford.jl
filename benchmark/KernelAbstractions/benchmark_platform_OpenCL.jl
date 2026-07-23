
#=============================================================================#
include("implementation/benchmark_platform.jl")

import pocl_jll
using OpenCL: CLArray, cl.devices, cl.platforms, cl.finish, cl.queue
const AT = CLArray
const path = "benchmarks/QuantumCliffordKAExt/OpenCL"

const can_run = any(length(devices(platform)) > 0 for platform in platforms())

if can_run
    synchronize() = finish(queue())
    benchmark_platform(synchronize, AT, path)
else
    @info "Unable to run OpenCL benchmark. No suitable device was found."
end
#=============================================================================#
