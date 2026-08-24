# GPU Extension Tests

This directory contains the tests for KernelAbstractions-supported hardware
accelerators. CUDA, ROCm, and OpenCL have separate test packages. Run one
backend at a time with `juliati` from the repository root,
for example `juliati . --filter ':cuda in tags'`. Buildkite assigns each tag to
a compatible runner.

# Noteworthy Details

- Synchronisation barriers are necessary before a test block because all
  accelerated operations are *asynchronous*.
- Some files contain workarounds for differences between this extension and the
  base package.
