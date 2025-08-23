# Adapt Extension

This directory contains the implementation for the Adapt derived functionality that enables smooth interoperation between objects residing in either the host or the device memory spaces. Please consult the JuliaGPU documentation for comprehensive information on how to setup and configure any specific device.

# Requirements

The following packages must be imported in order to activate this extension:
- [Adapt](https://github.com/JuliaGPU/Adapt.jl)
- [GPUArraysCore](https://github.com/JuliaGPU/GPUArrays.jl)

To benefit from the specialised execution path(s) available via dispatching to the underlying hardware accelerators, the [KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl) extension must be activated to enable those features.

# Noteworthy Details

- In order to support the utilisation of this extension in conjunction with the KernelAbstractions invocations, the bitwidth of the phase variable(s) must be compatible with the usage of atomic intrinsics. This is handled automatically by the adapt function calls but it introduces an incredibly minute discrepancy in the storage requirements depending on where the objects' memory is located.

# Warnings

The features provided herein remain an early and incomplete work-in-progress that is subject to continuous development. Bugs, missing features, and breaking changes are to be expected until such a time as when it is deemed suitable for official release. Consider this to be a thorough warning that **HERE BE DRAGONS**.
