# GPU Extension Tests

This directory contains the tests for KernelAbstractions-supported hardware
accelerators. CUDA, ROCm, and OpenCL have separate test environments. Resolve
and run one backend at a time from the repository root, for example:

```sh
project=test/gpu/cuda
julia --project="$project" -e 'import Pkg; Pkg.instantiate()'
juliati . --filter ':cuda in tags' \
  --env "JULIA_LOAD_PATH=@:$PWD/$project:@stdlib"
```

Buildkite assigns each backend to a compatible runner.

# Noteworthy Details

- Synchronisation barriers are necessary before a test block because all
  accelerated operations are *asynchronous*.
- Some files contain workarounds for differences between this extension and the
  base package.
