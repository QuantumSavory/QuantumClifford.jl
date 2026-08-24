# Test layout

`runtests.jl` is a compatibility entry point for `Pkg.test`, PkgEval, and the
standard GitHub test job. It runs only the test items stored directly in this
directory.

The specialized test suites are separate packages below this directory. Each
package has its own `Project.toml` and committed `Manifest.toml`, and depends on
the repository checkout of QuantumClifford and QECCore. The suites cover Aqua,
doctests, ECC, JET, Oscar, Python decoders, slow tests, and each GPU backend.

Install and run the TestItem application from the repository root:

```julia
import Pkg
Pkg.Apps.add("TestItemApp")
```

```sh
juliati .
```

Use `--filter` to select a package, test name, or tag. For example:

```sh
juliati . --filter 'package_name == "QuantumCliffordECCTests"'
juliati . --filter ':cuda in tags'
```

After changing dependencies in a specialized project, resolve its manifest
with that project active.
