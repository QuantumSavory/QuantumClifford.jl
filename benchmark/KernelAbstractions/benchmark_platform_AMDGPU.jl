include("benchmark_platform.jl")

using AMDGPU: ROCArray, synchronize, devices
AT = ROCArray
platform_name = "AMDGPU"

can_run = length(devices()) > 0

if can_run
	benchmark_platform(AT, synchronize; platform_name = platform_name)
else
	@info "Unable to benchmark $platform_name. No suitable device was found."
end
