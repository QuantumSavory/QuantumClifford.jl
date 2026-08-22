
#=============================================================================#
include("implementation/benchmark_platform.jl")

using CUDA: CuArray, devices, synchronize
const AT = CuArray
const path = "benchmarks/QuantumCliffordKAExt/CUDA"

const can_run = length(devices()) > 0

if can_run
    benchmark_platform(synchronize, AT, path)
else
    @info "Unable to run CUDA benchmark. No suitable device was found."
end
#=============================================================================#
