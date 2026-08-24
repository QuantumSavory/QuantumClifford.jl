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

```sh
juliati .
```

Use `--filter` to select a test name, file, or tag. For example:

```sh
juliati . --filter ':ecc in tags'
juliati . --filter ':cuda in tags'
```

Declare dependency and compatibility changes in the specialized project's
`Project.toml`.
