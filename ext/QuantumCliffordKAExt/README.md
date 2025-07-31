# KernelAbstractions Extension

This directory contains the implementation for the KernelAbstractions derived functionality that may be utilised in conjunction with any of the supported hardware accelerators. Please consult the JuliaGPU documentation for comprehensive information on how to setup and configure any specific device.

# Requirements

The following packages must be imported in order to activate this extension:
- [Atomix](https://github.com/JuliaConcurrent/Atomix.jl)
- [GPUArraysCore](https://github.com/JuliaGPU/GPUArrays.jl)
- [KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl)

In order to actually invoke its features, it is also pivotal to import the pertinent backend and ensure that the QuantumClifford defined objects are attached to storage residing within device memory. For convenience, an [Adapt](https://github.com/JuliaGPU/Adapt.jl) extension is also provided for handling the necessary conversion.

# Noteworthy Details

- It cannot be stressed enough that **ALL** the accelerated functionality is strictly *asynchronous* and that synchronisation barriers should be inserted as required. Please consult the relevant backend documentation for detailed instructions on this matter.
- Certain function calls are presently *synchronous* due to external limitations imposed by the toolchain dependencies. They should still be treated as *asynchronous* from a concurrency perspective as this will hopefully be resolved in a future release.
- Hardware limitations impose certain restrictions that are not present in the base package. Namely, the bitwidth of the phase variable(s) must be compatible with the usage of atomic intrinsics. The Adapt extension automatically handles this conversion but explicitly initialised objects must ensure their own compatibility.
- Wherever feasible, tuning parameters are exposed via `::Val` keyword arguments. Whilst the chosen defaults strive to be as performant as possible whilst maintaining generality, it may prove beneficial to further refine them to be more optimal for the underlying hardware.

# Warnings

The features provided herein remain an early and incomplete work-in-progress that is subject to continuous development. Bugs, missing features, and breaking changes are to be expected until such a time as when it is deemed suitable for official release. Consider this to be a thorough warning that **HERE BE DRAGONS**.
