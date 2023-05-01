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
            <a href="https://github.com/QuantumSavory/QuantumClifford.jl/actions?query=workflow%3ACI+branch%3Amaster"><img src="https://img.shields.io/github/actions/workflow/status/QuantumSavory/QuantumClifford.jl/ci.yml?branch=master" alt="GitHub Workflow Status"></a>
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

The only other simulator of similar performance I know of is [Stim](https://github.com/quantumlib/Stim). In particular, Stim implements convenient tracking of Pauli frames, that makes simulating the performance of error correcting codes blazingly fast (which are possible in QuantumClifford.jl, but no convenient interface is provided for that yet).

The "low level" functionality is of similar performance in Stim and QuantumClifford but different tradeoffs are made at the higher levels: to multiply in-place 1M-qubit Pauli operators Stim needs 16μs while QuantumClifford.jl needs 14μs. The difference is inconsequential and depends on compilers and hardware.

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

#### Gate application (500 CNOT gates on 1000 qubits) in 5 ms

```jldoctest
julia> @benchmark apply!(s, gate) setup=(s=random_stabilizer(1000); gate=tensor_pow(tCNOT,500))
BenchmarkTools.Trial: 931 samples with 1 evaluation.
 Range (min … max):  4.902 ms …   9.070 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     5.097 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   5.172 ms ± 319.591 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▃▇█▄▁▂▁▃▁                                                  
  ▃▇█████████▆▆▃▆▄▄▃▃▂▃▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂▁▂▂▁▂▂▁▁▁▁▁▁▁▃▁▁▂▁▁▁▁▂ ▃
  4.9 ms          Histogram: frequency by time         6.6 ms <

 Memory estimate: 6.89 KiB, allocs estimate: 49.
```

#### Sparse gate application to only specified qubits in a 1000 qubit tableau in 3 μs

```jldoctest
julia> @benchmark apply!(s, sCNOT(32,504)) setup=(s=random_stabilizer(1000))
BenchmarkTools.Trial: 10000 samples with 9 evaluations.
 Range (min … max):  2.602 μs …  12.860 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.934 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   3.048 μs ± 595.358 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

      ▂▆█▅                                                     
  ▁▁▂▅█████▅▃▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  2.6 μs          Histogram: frequency by time        5.53 μs <

 Memory estimate: 112 bytes, allocs estimate: 2.
```

#### Measuring a dense 1000 qubit Pauli operator in 70 μs

```jldoctest
julia> s=random_destabilizer(1000); p=random_pauli(1000);

julia> @benchmark project!(_s,_p) setup=(_s=copy(s);_p=copy(p)) evals=1
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  52.930 μs … 112.682 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     67.628 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   67.100 μs ±   5.624 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▂▁      ▃▄       ▆▆▄▁▁▂▃▃▃▄█▇▃▂▂▂▁▁▃▃▁         ▁          ▁▂ ▂
  ██▆▅▅▇▆▄██▇▅▅▄▄▄▆██████████████████████▇▇▇▇▇▆▇███▇▇█▇▇▇▆▇▇██ █
  52.9 μs       Histogram: log(frequency) by time        85 μs <

 Memory estimate: 480 bytes, allocs estimate: 4.
```

#### Measuring a single qubit in a 1000 qubit tableau in 50 μs

```jldoctest
julia> s=MixedDestabilizer(random_destabilizer(1000));

julia> @benchmark projectY!(_s,42) setup=(_s=copy(s)) evals=1
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  40.356 μs … 151.946 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     49.203 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   49.992 μs ±   6.442 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▃▁▄▄        ▄▇▅██▄▄▃▃▄▂           ▁ ▁▂▁                      ▂
  ████▇▅▅▄▅▄▅▅█████████████▇▇█▇▇▇████████▆▇▅▆▆▆▆▅▄▄▅▄▄▅▅▄▆▆▄▅▆ █
  40.4 μs       Histogram: log(frequency) by time      73.8 μs <

 Memory estimate: 432 bytes, allocs estimate: 3.
```

Benchmarks executed on a Ryzen Zen1 8-core CPU.

More detailed benchmarks can be seen at [github.com/QuantumSavory/QuantumCliffordBenchmarksLog](https://github.com/QuantumSavory/QuantumCliffordBenchmarksLog).
</details>
