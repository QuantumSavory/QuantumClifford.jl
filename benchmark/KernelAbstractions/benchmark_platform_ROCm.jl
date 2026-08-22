
#=============================================================================#
include("implementation/benchmark_platform.jl")

using AMDGPU: ROCArray, devices, synchronize
const AT = ROCArray
const path = "benchmarks/QuantumCliffordKAExt/ROCm"

const can_run = length(devices()) > 0

if can_run
    benchmark_platform(synchronize, AT, path)
else
    @info "Unable to run ROCm benchmark. No suitable device was found."
end
#=============================================================================#
