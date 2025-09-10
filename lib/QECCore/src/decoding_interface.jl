abstract type AbstractDecodingProblem end
abstract type AbstractNoiseModel end
abstract type AbstractDecodingScenario end

"""
    decoding_problem(code::AbstractCode, noise_model::AbstractNoiseModel, decoding_scenario::AbstractDecodingScenario)

Generate a decoding problem from a code, a noise model, and a decoding scenario.

### Inputs
- `code::AbstractCode`: The code to decode.
- `noise_model::AbstractNoiseModel`: The noise model to use.
- `decoding_scenario::AbstractDecodingScenario`: The decoding scenario to use.

### Outputs
- `problem::AbstractDecodingProblem`: The decoding problem.
"""
function decoding_problem end

abstract type AbstractDecodingSamples end
abstract type AbstractSampler end

"""
    sample(problem::AbstractDecodingProblem, sampler::AbstractSampler)

Sample from the decoding problem.

### Inputs
- `problem::AbstractDecodingProblem`: The decoding problem to sample from.
- `sampler::AbstractSampler`: The sampler to use.

### Outputs
- `samples::AbstractDecodingSamples`: The sampled syndrome.
"""
function sample end

abstract type AbstractSyndrome end
"""
    syndrome(samples::AbstractDecodingSamples)

Extract the syndrome from the samples.

### Inputs
- `samples::AbstractDecodingSamples`: The sampled data.

### Outputs
- `syndrome::AbstractSyndrome`: The extracted syndrome.
"""
function syndrome end

abstract type AbstractDecoder end

"""
    decode(problem::AbstractDecodingProblem, syndrome::AbstractSyndrome, decoder::AbstractDecoder)

Decode the syndrome using the decoder.

### Inputs
- `problem::AbstractDecodingProblem`: The decoding problem to decode.
- `syndrome::AbstractSyndrome`: The syndrome to decode.
- `decoder::AbstractDecoder`: The decoder to use.

### Outputs
- `decoding_result::AbstractMatrix`: The decoded result.
"""
function decode end

"""
    decoding_error_rate(problem::AbstractDecodingProblem, samples::AbstractDecodingSamples, decoding_result::AbstractMatrix)

Calculate the error rate of the decoding result.

### Inputs
- `problem::AbstractDecodingProblem`: The decoding problem to validate.
- `syndrome::AbstractSyndrome`: The syndrome to validate.
- `decoding_result::AbstractMatrix`: The decoded result to validate.

### Outputs
- `rate::Float64`: The error rate of the decoding result.
"""
function decoding_error_rate end
