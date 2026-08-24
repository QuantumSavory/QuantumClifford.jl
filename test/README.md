# QuantumClifford Tests

The root test suite uses ParallelTestRunner 2.8.1. Run commands from the package root. Positional arguments select tests by prefix; if no argument is given, the runner selects `general`.

The standard selectors are:

- `general`, including QEC, Oscar, and Python-decoder tests
- `slow`
- `jet`
- `KernelAbstractions/CUDA`
- `KernelAbstractions/ROCm`
- `KernelAbstractions/OpenCL`

Prepare the base environment before a local run:

```sh
julia --project=test -e 'using Pkg; Pkg.instantiate()'
```

The special test groups use standalone projects. Prepare only the projects required for the selected tests:

```sh
julia --project=test/projects/oscar -e 'using Pkg; Pkg.instantiate()'
julia --project=test/projects/python_decoders -e 'using Pkg; Pkg.instantiate()'
julia --project=test/projects/python_decoders -e 'using CondaPkg; CondaPkg.resolve()'
julia --project=test/projects/jet -e 'using Pkg; Pkg.instantiate()'
julia --project=test/projects/cuda -e 'using Pkg; Pkg.instantiate()'
julia --project=test/projects/rocm -e 'using Pkg; Pkg.instantiate()'
julia --project=test/projects/opencl -e 'using Pkg; Pkg.instantiate()'
```

Run a complete CPU group or a narrower prefix:

```sh
julia --project=test test/runtests.jl general
julia --project=test test/runtests.jl slow
julia --project=test test/runtests.jl jet
julia --project=test test/runtests.jl general/ecc
```

GPU groups must run serially because each file owns its backend worker:

```sh
julia --project=test test/runtests.jl --jobs=1 KernelAbstractions/CUDA
julia --project=test test/runtests.jl --jobs=1 KernelAbstractions/ROCm
julia --project=test test/runtests.jl --jobs=1 KernelAbstractions/OpenCL
```

Use `--list` to inspect the discovered tests without running them:

```sh
julia --project=test test/runtests.jl --list
```
