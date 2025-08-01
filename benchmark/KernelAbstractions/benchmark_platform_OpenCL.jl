include("implementation/benchmark_platform.jl")

using Dates: value, now, UNIXEPOCH
import pocl_jll
using OpenCL: CLArray, cl.devices, cl.platforms, cl.finish, cl.queue
const AT = CLArray
const path = "OpenCL_benchmark_" * string(value(now()) - UNIXEPOCH)

const can_run = any(length(devices(platform)) > 0 for platform in platforms())

if can_run
    synchronize() = finish(queue())
    benchmark_platform(AT, synchronize, path)
else
    @info "Unable to run OpenCL benchmark. No suitable device was found."
end
