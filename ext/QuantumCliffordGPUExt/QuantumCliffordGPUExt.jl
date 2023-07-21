module QuantumCliffordGPUExt

using CUDA
using QuantumClifford
import QuantumClifford: to_cpu, to_gpu, _apply!, apply!, pftrajectories

include("types.jl")
include("adapters.jl")
include("apply.jl")
include("pauli_frames.jl")
include("fastmemlayout.jl")

# todo hide _apply function later. only for internal use
end
