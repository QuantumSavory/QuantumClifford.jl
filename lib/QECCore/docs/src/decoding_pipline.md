# Quantum Error Correction Decoding Interface

## Overview

`decoding_interface.jl` defines the core abstract interfaces for quantum error correction decoding systems, providing a modular framework for handling decoding problems of quantum codes.

## Core Abstract Types

### Problem Definition
- `AbstractDecodingProblem`: Abstract base class for decoding problems
- `AbstractNoiseModel`: Abstract base class for noise models  
- `AbstractDecodingScenario`: Abstract base class for decoding scenarios

### Sampling and Data
- `AbstractDecodingSamples`: Abstract base class for decoding samples
- `AbstractSampler`: Abstract base class for samplers
- `AbstractSyndrome`: Abstract base class for syndromes

### Decoder
- `AbstractDecoder`: Abstract base class for decoders

## Main Interface Functions

### 1. Problem Generation
```julia
decoding_problem(code::AbstractCode, noise_model::AbstractNoiseModel, decoding_scenario::AbstractDecodingScenario) -> AbstractDecodingProblem
```
Generate a decoding problem from a quantum code, noise model, and decoding scenario.

### 2. Sampling
```julia
sample(problem::AbstractDecodingProblem, sampler::AbstractSampler) -> AbstractDecodingSamples
```
Sample from the decoding problem to generate decoding samples.

### 3. Syndrome Extraction
```julia
syndrome(samples::AbstractDecodingSamples) -> AbstractSyndrome
```
Extract syndrome from sampled data.

### 4. Decoding
```julia
decode(problem::AbstractDecodingProblem, syndrome::AbstractSyndrome, decoder::AbstractDecoder) -> AbstractMatrix
```
Decode the syndrome using the decoder and return the decoding result.

### 5. Error Rate Calculation
```julia
decoding_error_rate(problem::AbstractDecodingProblem, samples::AbstractDecodingSamples, decoding_result::AbstractMatrix) -> Float64
```
Calculate the error rate of the decoding result for performance evaluation.
