# KernelAbstractions Extension Tests

This directory contains the test suite for the range of KernelAbstractions-supported hardware accelerators. These tests require compatible hardware and use separate projects for CUDA, ROCm, and OpenCL. Select one backend by its positional test prefix and use one ParallelTestRunner worker. For example:

```sh
julia --project=test test/runtests.jl --jobs=1 KernelAbstractions/CUDA
```

# Noteworthy Details

- Synchronisation barriers are necessary prior to any test block in order to ensure their validity, given that all accelerated functionality is *asynchronous*.
- Some files require additional massaging and roundabout workarounds in order to perform their tests, which is due to a certain level of feature disparity between this extension and the base package. This will hopefully be resolved in a future release.
