module QuantumCliffordGPUExt

using CUDA
using QuantumClifford
import QuantumClifford: to_cpu, to_gpu, _apply!, apply!, pftrajectories, fastcolumn, fastrow, applynoise!

include("utils.jl")
include("types.jl")
include("adapters.jl")
include("apply.jl")
include("pauli_frames.jl")
include("fastmemlayout.jl")
include("apply_noise.jl")

# todo hide _apply function later. only for internal use
end
