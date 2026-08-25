# Test layout

`runtests.jl` is a compatibility entry point for `Pkg.test`, PkgEval, and the
standard GitHub test job. It runs only the test items stored directly in this
directory.

The specialized test suites use separate environments below this directory.
Each environment has its own `Project.toml` and depends on the repository
checkout of QuantumClifford and QECCore. Dependencies resolve on each run so
that these suites test current compatible package versions. The suites cover
Aqua, doctests, ECC, JET, Oscar, Python decoders, slow tests, and each GPU
backend.

Install and run the TestItem application from the repository root:

```julia
import Pkg
Pkg.Apps.add("TestItemApp")
```

Run the base items directly:

```sh
juliati . --filter 'dirname(filename) == realpath("test")'
```

Resolve a specialized environment and expose it only to that suite's test
processes. For example:

```sh
project=test/ecc
julia --project="$project" -e 'import Pkg; Pkg.instantiate()'
juliati . \
  --filter 'startswith(filename, string(realpath("test/ecc"), Base.Filesystem.path_separator))' \
  --env "JULIA_LOAD_PATH=@:$PWD/$project:@stdlib"
```

Buildkite runs each non-hardware environment as a separate matrix job and each
GPU environment on its hardware runner. Runtime manifests remain uncommitted,
so every job resolves current compatible dependency versions. Declare
dependency and compatibility changes in the corresponding `Project.toml`.
