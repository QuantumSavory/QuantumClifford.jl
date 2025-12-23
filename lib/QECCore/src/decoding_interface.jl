abstract type AbstractDecodingProblem end
abstract type AbstractDecodingScenario end

"""
    decoding_problem(code::AbstractCode, decoding_scenario::AbstractDecodingScenario)

Generate a decoding problem from a code and a decoding scenario.

### Inputs
- `code::AbstractCode`: The code to decode.
- `decoding_scenario::AbstractDecodingScenario`: The decoding scenario to use.

### Outputs
- `problem::AbstractDecodingProblem`: The decoding problem.
"""
function decoding_problem end

abstract type AbstractDecodingSamples end
abstract type AbstractSampler end

"""
    sample(problem::AbstractDecodingProblem, num_samples::Int)
    sample(problem::AbstractDecodingProblem, num_samples::Int, sampler::AbstractSampler)

Sample from the decoding problem. If the sampler is not provided, the default sampler will be used.

### Inputs
- `problem::AbstractDecodingProblem`: The decoding problem to sample from.
- `num_samples::Int`: The number of samples to generate.
- `sampler::AbstractSampler`: The sampler to use.

### Outputs
- `samples::AbstractDecodingSamples`: The sampled data.
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
abstract type AbstractDecodingResult end
"""
    decode(problem::AbstractDecodingProblem, syndrome::AbstractSyndrome, decoder::AbstractDecoder)

Decode the syndrome using the decoder.

### Inputs
- `problem::AbstractDecodingProblem`: The decoding problem to decode.
- `syndrome::AbstractSyndrome`: The syndrome to decode.
- `decoder::AbstractDecoder`: The decoder to use.

### Outputs
- `decoding_result::AbstractDecodingResult`: The decoded result.
"""
function decode end

"""
    decoding_error_rate(problem::AbstractDecodingProblem, samples::AbstractDecodingSamples, decoding_result::AbstractDecodingResult)

Calculate the error rate of the decoding result.

### Inputs
- `problem::AbstractDecodingProblem`: The decoding problem to validate.
- `samples::AbstractDecodingSamples`: The samples to validate.
- `decoding_result::AbstractDecodingResult`: The decoded result to validate.

### Outputs
- `rate::Float64`: The error rate of the decoding result.
"""
function decoding_error_rate end
