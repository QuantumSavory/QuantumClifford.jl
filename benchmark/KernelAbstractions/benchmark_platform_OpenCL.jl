include("benchmark_platform.jl")

import pocl_jll
using OpenCL: CLArray, cl.devices, cl.platforms, cl.finish, cl.queue
AT = CLArray
platform_name = "OpenCL"

can_run = any(length(devices(platform)) > 0 for platform in platforms())

if can_run
	synchronize() = finish(queue())
	benchmark_platform(AT, synchronize; platform_name = platform_name)
else
	@info "Unable to benchmark $platform_name. No suitable device was found."
end
