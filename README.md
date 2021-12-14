# QuantumClifford.jl

[![Documentation of latest stable version](https://img.shields.io/badge/docs-stable-blue.svg)](https://krastanov.github.io/QuantumClifford.jl/stable)
[![Documentation of dev version](https://img.shields.io/badge/docs-dev-blue.svg)](https://krastanov.github.io/QuantumClifford.jl/dev)
[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/Krastanov/QuantumClifford.jl/CI)](https://github.com/Krastanov/QuantumClifford.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![Test coverage from codecov](https://img.shields.io/codecov/c/gh/Krastanov/QuantumClifford.jl?label=codecov)](https://codecov.io/gh/Krastanov/QuantumClifford.jl)
<!--[![Build status](https://api.travis-ci.com/Krastanov/QuantumClifford.jl.svg?branch=master)](https://travis-ci.com/Krastanov/QuantumClifford.jl)-->
<!--[![Test coverage from coveralls](https://img.shields.io/coveralls/github/Krastanov/QuantumClifford.jl?label=coveralls)](https://coveralls.io/r/Krastanov/QuantumClifford.jl?branch=master)-->

A Julia package for working with quantum stabilizer states and Clifford circuits
that act on them. The package is already very fast for the majority of common operations, but there are still many low-hanging fruits performance-wise. See the detailed [suggested readings & references page](https://krastanov.github.io/QuantumClifford.jl/dev/references/#Suggested-reading) for background on the various algorithms.

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

The code is vectorized and multithreaded.


## Quick Benchmarks

Fast, in-place, allocation free implementations.

#### Comparison against other Clifford simulators

The only other simulator of similar performance I know of is [Stim](https://github.com/quantumlib/Stim). In particular, Stim implements convenient tracking of Pauli frames, that makes simulating the performance of error correcting codes blazingly fast (which are possible in QuantumClifford.jl, but no convenient interface is provided for that yet).

The "low level" functionality is of similar performance in Stim and QuantumClifford but different tradeoffs are made at the higher levels: to multiply in-place 1M-qubit Pauli operators Stim needs 16us while QuantumClifford.jl needs 20us.

Of note is that Stim achieved this performance through high-quality C++ SIMD code of significant sophistication, while QuantumClifford.jl is implemented in pure Julia.

#### Multiplying two 1 gigaqubit Paulis in 32 ms

```jldoctest
julia> a = random_pauli(1_000_000_000);
julia> b = random_pauli(1_000_000_000);
julia> @benchmark QuantumClifford.mul_left!(a,b)
BenchmarkTools.Trial: 155 samples with 1 evaluation.
 Range (min … max):  32.074 ms … 32.425 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     32.246 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   32.247 ms ± 63.427 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

                  ▃  ▃▃ ▄ ▆▄▄▄██▃ ▃▄▁▆█▃▃ ▃      ▁             
  ▄▁▁▄▁▁▄▆▁▁▄▆▄▆▇▇█▄▄██▄█▆███████▇███████▆█▆▄▄▄▁▄█▁▄▄▁▄▁▁▁▁▁▄ ▄
  32.1 ms         Histogram: frequency by time        32.4 ms <

 Memory estimate: 0 bytes, allocs estimate: 0.
```

#### Canonicalization of a random 1000-qubit stabilizer in 22 ms

```jldoctest
julia> @benchmark canonicalize!(s) setup=(s=random_stabilizer(1000))
BenchmarkTools.Trial: 226 samples with 1 evaluation.
 Range (min … max):  21.938 ms …  22.680 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     22.025 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   22.057 ms ± 115.247 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▂▂ █▃▃▂                                                     
  ▄▇███████▆▇▆█▆▄▄▄▄▄▅▄▃▃▁▄▃▃▃▃▃▃▁▁▁▁▁▃▁▁▁▁▁▃▁▃▁▁▁▁▃▁▁▁▁▁▁▁▁▃▃ ▃
  21.9 ms         Histogram: frequency by time         22.6 ms <

 Memory estimate: 32 bytes, allocs estimate: 1.
```

#### Gate application (500 CNOT gates on 1000 qubits) in 6 ms

```jldoctest
julia> @benchmark apply!(s, gate) setup=(s=random_stabilizer(1000); gate=tensor_pow(CNOT,500))
BenchmarkTools.Trial: 801 samples with 1 evaluation.
 Range (min … max):  5.758 ms …  10.479 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     5.800 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   6.107 ms ± 624.182 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █▆▃                               ▂                          
  ███▇▆▅▅▅▅▅▆▅▄▁▁▆▄▇▆▆▆▇▅▅▅▅▅▆▆▆▅▇███▇▆▇▇▅▅▁▅▇▅▅▁▅▅▁▁▁▁▁▄▁▁▁▅ ▇
  5.76 ms      Histogram: log(frequency) by time      8.22 ms <

 Memory estimate: 52.62 KiB, allocs estimate: 2592.
```

#### Sparse gate application to only specified qubits in a 1000 qubit tableau in 4 microseconds

```jldoctest
julia> @benchmark apply!(s, sCNOT(32,504)) setup=(s=random_stabilizer(1000))
BenchmarkTools.Trial: 10000 samples with 9 evaluations.
 Range (min … max):  3.373 μs … 252.630 μs  ┊ GC (min … max): 0.00% … 53.27%
 Time  (median):     3.766 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   3.892 μs ±   2.525 μs  ┊ GC (mean ± σ):  0.35% ±  0.53%

        ▃▆█▅▁                                                  
  ▁▁▁▂▃▇█████▅▃▂▂▃▃▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  3.37 μs         Histogram: frequency by time        6.07 μs <

 Memory estimate: 96 bytes, allocs estimate: 2.
```

Benchmarks executed on a Ryzen Zen1 8-core CPU.