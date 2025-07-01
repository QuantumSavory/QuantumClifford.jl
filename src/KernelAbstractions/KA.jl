
#=============================================================================#
using GPUArraysCore: AbstractGPUArray
import KernelAbstractions as KA
# Resolves issue due to KA comparing against the literal Symbol("@Const").
using KernelAbstractions: @Const
import Atomix

# Most ideal type for shared reductions due to avoiding bank conflicts.
const DeviceUnsigned = UInt32
# Reasonable size that is generally ideal for most vendors and usecases.
const default_block_size = 256
# Ameliorate overhead and enhance performance by doing more work per thread.
const default_batch_size = 32

include("utils.jl")
include("mul_leftright.jl")
#=============================================================================#
