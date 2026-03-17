# KernelAbstractions Extension Tests

This directory contains the test suite for the range of KernelAbstractions-supported hardware accelerators. These tests are disabled by default since compatible hardware may not necessarily exist on all systems. Only a single platform may be tested in any individual run and it must be explicitly enabled by having the environment set `[CUDA, ROCm, OpenCL]_TEST=true` as desired.

# Noteworthy Details

- Synchronisation barriers are necessary prior to any test block in order to ensure their validity, given that all accelerated functionality is *asynchronous*.
- Some files require additional massaging and roundabout workarounds in order to perform their tests, which is due to a certain level of feature disparity between this extension and the base package. This will hopefully be resolved in a future release.
