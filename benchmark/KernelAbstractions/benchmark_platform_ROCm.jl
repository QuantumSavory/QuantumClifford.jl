include("implementation/benchmark_platform.jl")

using Dates: value, now, UNIXEPOCH
using AMDGPU: ROCArray, devices, synchronize
const AT = ROCArray
const path = "ROCm_benchmark_" * string(value(now()) - UNIXEPOCH)

const can_run = length(devices()) > 0

if can_run
    benchmark_platform(AT, synchronize, path)
else
    @info "Unable to run ROCm benchmark. No suitable device was found."
end
