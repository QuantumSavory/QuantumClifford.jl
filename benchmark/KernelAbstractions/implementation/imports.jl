
#=============================================================================#
# Required for QuantumCliffordKAExt.
import Atomix, GPUArraysCore, KernelAbstractions

# Required for QuantumCliffordAdaptExt.
using Adapt: adapt

using BenchmarkTools: @belapsed

using Dates: value, now, UNIXEPOCH, Second

# Assists in reducing resource demands.
using GPUArrays: AllocCache, @cached, unsafe_free!

# Utilised for extrapolating the runtime of exceedingly long benchmarks.
using LsqFit: curve_fit

using Plots: plot, savefig

using QuantumClifford
# This must be done explicitly as they are not exported.
using QuantumClifford: Tableau, AbstractStabilizer, random_tableau
#=============================================================================#
