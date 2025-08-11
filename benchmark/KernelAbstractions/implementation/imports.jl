# Required for QuantumCliffordKAExt.
import Atomix, GPUArraysCore, KernelAbstractions

using BenchmarkTools: @belapsed
using Dates: value, now, UNIXEPOCH
# Assists in reducing resource demands.
using GPUArrays: AllocCache, @cached, unsafe_free!
using Plots: plot, savefig
using QuantumClifford
