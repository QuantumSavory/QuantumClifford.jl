# QuantumClifford.jl

[![Documentation of latest stable version](https://img.shields.io/badge/docs-stable-blue.svg)](https://krastanov.github.io/QuantumClifford.jl/stable)
[![Documentation of dev version](https://img.shields.io/badge/docs-dev-blue.svg)](https://krastanov.github.io/QuantumClifford.jl/dev)
[![Build Status](https://github.com/Krastanov/QuantumClifford.jl/workflows/CI/badge.svg)](https://github.com/Krastanov/QuantumClifford.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![Build status](https://api.travis-ci.com/Krastanov/QuantumClifford.jl.svg?branch=master)](https://travis-ci.com/Krastanov/QuantumClifford.jl)
[![Test coverage from codecov](https://img.shields.io/codecov/c/gh/Krastanov/QuantumClifford.jl?label=codecov)](https://codecov.io/gh/Krastanov/QuantumClifford.jl)
[![Test coverage from coveralls](https://img.shields.io/coveralls/github/Krastanov/QuantumClifford.jl?label=coveralls)](https://coveralls.io/r/Krastanov/QuantumClifford.jl?branch=master)

A Julia package for working with quantum stabilizer states and Clifford circuits
that act on them.

The package is still in an alpha state. It is already very fast for the majority of common operations, but there are still many low-hanging fruits performance-wise.

To install it use:

```julia
] add QuantumClifford
```

Works efficiently with
[pure](https://krastanov.github.io/QuantumClifford/dev/manual/#Stabilizers-1) and
[mixed stabilizer](https://krastanov.github.io/QuantumClifford/dev/mixed/#Mixed-Stabilizer-States-1)
states of thousands of qubits
as well as
[sparse or dense Clifford operations](https://krastanov.github.io/QuantumClifford/dev/manual/#Clifford-Operators-1)
acting upon them.

Provides
[canonicalization](https://krastanov.github.io/QuantumClifford/dev/manual/#Canonicalization-of-Stabilizers-1),
[projection](https://krastanov.github.io/QuantumClifford/dev/manual/#Projective-Measurements-1), and
[generation](https://krastanov.github.io/QuantumClifford/dev/manual/#Generating-a-Pauli-Operator-with-Stabilizer-Generators-1) operations,
as well as
[partial traces](https://krastanov.github.io/QuantumClifford/dev/manual/#Partial-Traces-1).

```jldoctest
julia> P"X" * P"Z"
-iY

julia> P"X" ⊗ P"Z"
+ XZ

julia> S"-XX
         +ZZ"
- XX
+ ZZ

julia> CNOT * S"-XX
                +ZZ"
- X_
+ _Z
```


## Quick Benchmarks

Fast, in-place, (mostly) allocation free implementations.

#### Canonicalization of a random 100-qubit stabilizer

```jldoctest
julia> @benchmark canonicalize!(s) setup=(s=random_stabilizer(100))
BenchmarkTools.Trial:
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     139.927 μs (0.00% GC)
  median time:      161.090 μs (0.00% GC)
  mean time:        165.824 μs (0.00% GC)
  maximum time:     278.056 μs (0.00% GC)
  --------------
  samples:          376
  evals/sample:     1
```

#### Gate application (50 CNOT gates on 100 qubits)

```jldoctest
julia> @benchmark apply!(s, gate) setup=(s=random_stabilizer(100); gate=tensor_pow(CNOT,50))
BenchmarkTools.Trial:
  memory estimate:  95.09 KiB
  allocs estimate:  6045
  --------------
  minimum time:     416.188 μs (0.00% GC)
  median time:      482.586 μs (0.00% GC)
  mean time:        492.410 μs (0.00% GC)
  maximum time:     687.202 μs (0.00% GC)
  --------------
  samples:          101
  evals/sample:     1
```

#### Sparse gate application to only specified qubits

```jldoctest
julia> @benchmark apply!(s, CNOT, [32,54]) setup=(s=random_stabilizer(100))
BenchmarkTools.Trial:
  memory estimate:  2.67 KiB
  allocs estimate:  126
  --------------
  minimum time:     12.425 μs (0.00% GC)
  median time:      15.268 μs (0.00% GC)
  mean time:        15.713 μs (0.00% GC)
  maximum time:     32.185 μs (0.00% GC)
  --------------
  samples:          407
  evals/sample:     1
```
