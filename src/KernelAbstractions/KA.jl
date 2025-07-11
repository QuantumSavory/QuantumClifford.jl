
#=============================================================================#
using GPUArraysCore: AbstractGPUArray
import KernelAbstractions as KA
# Resolves issue due to KA comparing against the literal Symbol("@Const").
using KernelAbstractions: @Const
import Atomix

include("definitions.jl")
include("utilities.jl")
include("mul_leftright.jl")
#=============================================================================#
