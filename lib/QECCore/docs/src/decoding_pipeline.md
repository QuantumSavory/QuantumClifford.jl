# Quantum Error Correction Decoding Interface

## Overview

`decoding_interface.jl` defines the core abstract interfaces for quantum error correction decoding systems, providing a modular framework for handling decoding problems of quantum codes.

## Core Abstract Types

### Problem Definition
- `AbstractDecodingProblem`: Abstract base class for decoding problems
- `AbstractDecodingScenario`: Abstract base class for decoding scenarios

### Sampling and Data
- `AbstractDecodingSamples`: Abstract base class for decoding samples
- `AbstractSampler`: Abstract base class for samplers
- `AbstractSyndrome`: Abstract base class for syndromes

### Decoder
- `AbstractDecoder`: Abstract base class for decoders
- `AbstractDecodingResult`: Abstract base class for decoding results

## Main Interface Functions

### 1. Problem Generation
```julia
decoding_problem(code::AbstractCode, decoding_scenario::AbstractDecodingScenario) -> AbstractDecodingProblem
```
Generate a decoding problem from a quantum code and a decoding scenario.

### 2. Sampling
```julia
sample(problem::AbstractDecodingProblem) -> AbstractDecodingSamples
sample(problem::AbstractDecodingProblem, sampler::AbstractSampler) -> AbstractDecodingSamples
```
Sample from the decoding problem to generate decoding samples. If the sampler is not provided, the default sampler will be used.

### 3. Syndrome Extraction
```julia
syndrome(samples::AbstractDecodingSamples) -> AbstractSyndrome
```
Extract syndrome from sampled data.

### 4. Decoding
```julia
decode(problem::AbstractDecodingProblem, syndrome::AbstractSyndrome, decoder::AbstractDecoder) -> AbstractDecodingResult
```
Decode the syndrome using the decoder and return the decoding result.

### 5. Error Rate Calculation
```julia
decoding_error_rate(problem::AbstractDecodingProblem, samples::AbstractDecodingSamples, decoding_result::AbstractDecodingResult) -> Float64
```
Calculate the error rate of the decoding result for performance evaluation.

## Examples
The following examples demonstrate the typical workflow using the decoding interface. These are conceptual examples showing the general pattern of usage and are not be implemented yet.

### Code Capacity Noise
Code capacity noise is the simplest decoding scenario where noise only affects data qubits, and syndrome measurements are assumed to be perfect (no measurement errors). This scenario is useful for understanding the fundamental error correction capability of a code.

```julia
# Step 1: Create a code capacity decoding scenario
code = Steane7()
scenario = CodeCapacityScenario(error_rate = 0.1)
problem = decoding_problem(code, scenario)

# Step 2: Sample from the problem
samples = sample(problem)  # Uses default sampler

# Step 3: Extract syndrome
syndrome_data = syndrome(samples)

# Step 4: Decode using a decoder
decoder = TableDecoder(code)
result = decode(problem, syndrome_data, decoder)

# Step 5: Evaluate decoding performance
error_rate = decoding_error_rate(problem, samples, result)
# Checks if the correction leads to a logical error
```

In code capacity noise, the decoding problem is straightforward: given a perfect syndrome `s = H * e` (where `H` is the parity check matrix and `e` is the error pattern), find the most likely error pattern `e'` that satisfies `H * e' = s`.

### Circuit-level Decoding

Circuit-level decoding is a more realistic scenario where noise affects both data qubits and the syndrome measurement circuit. This includes measurement errors, gate errors, and other imperfections in the measurement process.

```julia
code = Toric(3,3)
scenario = CircuitLevelScenario(
    after_clifford_depolarization=0.01,
    after_reset_flip_probability=0.01,
    before_measure_flip_probability=0.01,
    before_round_data_depolarization=0.01,
)
problem = decoding_problem(code, scenario)
```
Since the struct of problem is also a `DetectorModelProblem`, the sampling and decoding process are the same as the code capacity noise.
