include("implementation/benchmark_platform.jl")

using CUDA: CuArray, devices, synchronize
const AT = CuArray
const path = "QuantumClifford_benchmarks/CUDA"

const can_run = length(devices()) > 0

if can_run
    benchmark_platform(AT, synchronize, path)
else
    @info "Unable to run CUDA benchmark. No suitable device was found."
end
