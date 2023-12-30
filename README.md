# QuantumClifford.jl

<table>
    <tr>
        <td>Documentation</td>
        <td>
            <a href="https://quantumsavory.github.io/QuantumClifford.jl/stable"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Documentation of latest stable version"></a>
            <a href="https://quantumsavory.github.io/QuantumClifford.jl/dev"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation of dev version"></a>
        </td>
    </tr><tr></tr>
    <tr>
        <td>Continuous integration</td>
        <td>
            <a href="https://github.com/QuantumSavory/QuantumClifford.jl/actions?query=workflow%3ACI+branch%3Amaster"><img src="https://github.com/QuantumSavory/QuantumClifford.jl/actions/workflows/ci.yml/badge.svg" alt="GitHub Workflow Status"></a>
            <a href="https://buildkite.com/quantumsavory/quantumclifford"><img src="https://badge.buildkite.com/8ef137151415f29c03544c5b7963f6bc6afc1022f29cfc072a.svg?branch=master" alt="Buildkite Workflow Status"></a>
        </td>
    </tr><tr></tr>
    <tr>
        <td>Code coverage</td>
        <td>
            <a href="https://codecov.io/gh/QuantumSavory/QuantumClifford.jl"><img src="https://img.shields.io/codecov/c/gh/QuantumSavory/QuantumClifford.jl?label=codecov" alt="Test coverage from codecov"></a>
        </td>
    </tr><tr></tr>
    <tr>
        <td>Static analysis with</td>
        <td>
            <a href="https://github.com/aviatesk/JET.jl"><img src="https://img.shields.io/badge/JET.jl-%E2%9C%88%EF%B8%8F-9cf" alt="JET static analysis"></a>
            <a href="https://github.com/JuliaTesting/Aqua.jl"><img src="https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg" alt="Aqua QA"></a>
        </td>
    </tr>
</table>


A Julia package for working with quantum stabilizer states and Clifford circuits
that act on them. Graphs states are also supported. The package is already very fast for the majority of common operations, but there are still many low-hanging fruits performance-wise. See the detailed [suggested readings & references page](https://quantumsavory.github.io/QuantumClifford.jl/dev/references/#Suggested-reading) for background on the various algorithms.

To install it use:

```julia
] add QuantumClifford
```

Works efficiently with
[pure](https://quantumsavory.github.io/QuantumClifford.jl/dev/manual/#Stabilizers-1) and
[mixed stabilizer](https://quantumsavory.github.io/QuantumClifford.jl/dev/mixed/#Mixed-Stabilizer-States-1)
states of thousands of qubits
as well as
[sparse or dense Clifford operations](https://quantumsavory.github.io/QuantumClifford.jl/dev/manual/#Clifford-Operators-1)
acting upon them.

Implements [Pauli frames](https://quantumsavory.github.io/QuantumClifford.jl/dev/ecc_example_sim/) for fast sampling.

Provides
[canonicalization](https://quantumsavory.github.io/QuantumClifford.jl/dev/manual/#Canonicalization-of-Stabilizers-1),
[projection](https://quantumsavory.github.io/QuantumClifford.jl/dev/manual/#Projective-Measurements-1), and
[generation](https://quantumsavory.github.io/QuantumClifford.jl/dev/manual/#Generating-a-Pauli-Operator-with-Stabilizer-Generators-1) operations,
as well as
[partial traces](https://quantumsavory.github.io/QuantumClifford.jl/dev/manual/#Partial-Traces-1).

```jldoctest
julia> P"X" * P"Z"
-iY

julia> P"X" ⊗ P"Z"
+ XZ

julia> S"-XX
         +ZZ"
- XX
+ ZZ

julia> tCNOT * S"-XX
                 +ZZ"
- X_
+ _Z
```

The code is vectorized and multithreaded.

Fast, in-place, allocation free implementations.

<details>
    <summary>Quick Benchmarks (click to expand)</summary>

#### Comparison against other Clifford simulators

The only other simulator of similar performance I know of is [Stim](https://github.com/quantumlib/Stim).

The "low level" functionality is of similar performance in Stim and QuantumClifford but different tradeoffs are made at the higher levels: to multiply in-place 1M-qubit Pauli operators Stim and QuantumClifford.jl both need around 15μs. The difference is inconsequential and depends on compilers and hardware.

Of note is that Stim achieved this performance through high-quality C++ SIMD code of significant sophistication, while QuantumClifford.jl is implemented in pure and simple Julia.

#### Multiplying two 1 gigaqubit Paulis in 13 ms

```jldoctest
julia> a = random_pauli(1_000_000_000);
julia> b = random_pauli(1_000_000_000);
julia> @benchmark QuantumClifford.mul_left!(a,b)
BenchmarkTools.Trial: 373 samples with 1 evaluation.
 Range (min … max):  13.209 ms …  14.304 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     13.355 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   13.427 ms ± 173.503 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 0 bytes, allocs estimate: 0.
```

#### Canonicalization of a random 1000-qubit stabilizer in 9 ms

```jldoctest
julia> @benchmark canonicalize!(s) setup=(s=random_stabilizer(1000))
BenchmarkTools.Trial: 6 samples with 1 evaluation.
 Range (min … max):  8.516 ms …  8.614 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     8.536 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   8.550 ms ± 35.883 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 0 bytes, allocs estimate: 0.
```

#### Dense tableaux multiplication (tensor product of 500 CNOT gates acting 1000 qubits) in 17 ms

```jldoctest
julia> @benchmark apply!(s, gate) setup=(s=random_stabilizer(1000); gate=tensor_pow(tCNOT,500))
BenchmarkTools.Trial: 6 samples with 1 evaluation.
 Range (min … max):  16.879 ms … 17.064 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     17.010 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   16.997 ms ± 63.050 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 800 bytes, allocs estimate: 4.
```

#### Sparse gate application to only specified qubits in a 1000 qubit tableau in 3 μs

```jldoctest
julia> @benchmark apply!(s, sCNOT(32,504)) setup=(s=random_stabilizer(1000))
BenchmarkTools.Trial: 6 samples with 8 evaluations.
 Range (min … max):  2.867 μs …   3.228 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     3.043 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   3.049 μs ± 119.106 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 0 bytes, allocs estimate: 0.
```

#### Measuring a dense 1000 qubit Pauli operator in 18 μs

```jldoctest
julia> s=random_destabilizer(1000); p=random_pauli(1000);

julia> @benchmark project!(_s,_p) setup=(_s=copy(s);_p=copy(p)) evals=1
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  17.753 μs … 39.444 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     21.971 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   21.893 μs ±  2.234 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 480 bytes, allocs estimate: 4.
```

#### Measuring a single qubit in a 1000 qubit tableau in 15 μs

```jldoctest
julia> s=MixedDestabilizer(random_destabilizer(1000));

julia> @benchmark projectY!(_s,42) setup=(_s=copy(s)) evals=1
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  15.379 μs … 37.630 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     16.912 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   17.120 μs ±  1.335 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

 Memory estimate: 464 bytes, allocs estimate: 5.
```

Benchmarks executed on a single thread on Ryzen Zen4 16-core CPU:

```
julia> versioninfo()
Julia Version 1.9.1
Commit 147bdf428cd (2023-06-07 08:27 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 32 × AMD Ryzen 9 7950X 16-Core Processor
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-14.0.6 (ORCJIT, znver3)
  Threads: 1 on 32 virtual cores
```

More detailed benchmarks can be seen at [github.com/QuantumSavory/QuantumCliffordBenchmarksLog](https://github.com/QuantumSavory/QuantumCliffordBenchmarksLog).
</details>
