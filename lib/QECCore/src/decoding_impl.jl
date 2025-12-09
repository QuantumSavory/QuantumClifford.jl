abstract type AbstractBitNoiseModel end
"""
    $TYPEDEF

Represents a probabilistic error model for a system of bits using a factored 
distribution approach. The model decomposes the joint probability distribution 
over all bits into smaller, manageable distributions over subsets of bits.

### Fields
    $TYPEDFIELDS
"""
struct FactoredBitNoiseModel{VT<:AbstractArray{Float64}} <: AbstractBitNoiseModel
    """The number of bits in the noise model."""
    num_bits::Int
    """Dictionary defining the factored probability distributions.

    - **Keys**: `Vector{Int}` specifying which subset of bits the distribution covers
      (e.g., `[1,3]` means this distribution covers bits 1 and 3)
    - **Values**: `VT` multi-dimensional probability array where:
        - Dimensionality = length of the key vector
        - Each dimension has size 2 (representing classical bit states 0/1 or error/no-error)
        - Values are non-negative probabilities summing to 1

    Example: `[1,2] => [0.4 0.3; 0.2 0.1]` represents:
        - P(bit1=0, bit2=0) = 0.4
        - P(bit1=0, bit2=1) = 0.3
        - P(bit1=1, bit2=0) = 0.2
        - P(bit1=1, bit2=1) = 0.1"""
    probabilities::Dict{Vector{Int},VT}
    function FactoredBitNoiseModel(num_bits::Int, probabilities::Dict{Vector{Int},VT}) where VT<:AbstractArray{Float64}
        for (error_bits, prob_array) in probabilities
            @assert maximum(error_bits) <= num_bits "Maximum element in error bits $error_bits is $(maximum(error_bits)), but num_bits is $num_bits"
            expected_size = (fill(2, length(error_bits))...,)
            actual_size = size(prob_array)
            @assert size(prob_array) == expected_size "The dimension of the probability array must match the length of the error bits vector, and each dimension must be size 2. Expected: $expected_size, Got: $actual_size"
            @assert sum(prob_array) == 1 "The sum of the probabilities must be 1, but got $(sum(prob_array))"
        end
        missing_variables = setdiff(1:num_bits, unique(vcat(keys(probabilities)...)))
        @assert isempty(missing_variables) "The following variables are not used in the error model: $missing_variables"
        new{VT}(num_bits, probabilities)
    end
end

"""
    $TYPEDSIGNATURES

Create an error model from a vector of error probabilities. The `i`-th element of the vector is the error probability of bit `i`.

See also: [`depolarization_error_model`](@ref)
"""
function FactoredBitNoiseModel(probabilities::Vector{Float64})
    num_bits = length(probabilities)
    probabilities = Dict([i] => [1.0 - probabilities[i], probabilities[i]] for i in 1:num_bits)
    return FactoredBitNoiseModel(num_bits, probabilities)
end

"""
    $TYPEDSIGNATURES

Create a depolarization error model from a vector of error probabilities. 
The `i`-th element of the vector is the depolarization probability of qubit `i`.
"""
function depolarization_error_model(pvec::Vector{Float64}, qubit_num::Int)
    @assert length(pvec) == qubit_num "The length of the vector of error probabilities must match the number of qubits"
    probabilities = Dict([i, i + qubit_num] => [1.0-pvec[i] pvec[i]/3; pvec[i]/3 pvec[i]/3] for i in 1:qubit_num)
    return FactoredBitNoiseModel(2 * qubit_num, probabilities)
end

"""
    $TYPEDSIGNATURES

Create a depolarization error model from a depolarization probability. All qubits have the same depolarization probability.
"""
depolarization_error_model(p::Float64, qubit_num::Int) = depolarization_error_model(fill(p, qubit_num), qubit_num)

"""
    $TYPEDSIGNATURES

Check if the error model is a vector error model.
"""
isvector(em::FactoredBitNoiseModel) = all(length(key) == 1 for key in keys(em.probabilities))

"""
    $TYPEDSIGNATURES

Check if the error model is an independent error model. Independent error model means that there is no overlap between the error bits of different probability distributions.
"""
function isindependent(em::FactoredBitNoiseModel)
    keys_list = collect(keys(em.probabilities))
    for i in 1:length(keys_list)-1
        for j in i+1:length(keys_list)
            if !isempty(intersect(keys_list[i], keys_list[j]))
                return false
            end
        end
    end
    return true
end

"""
    $TYPEDEF

Represents an independent error model for a system of bits.
"""
struct IndependentBitNoiseModel <: AbstractBitNoiseModel
    """The error probability of the bits"""
    probabilities::Vector{Float64}
    function IndependentBitNoiseModel(probabilities::Vector{Float64})
        @assert all(0 .<= probabilities .<= 1) "The error probabilities must be between 0 and 1"
        new(probabilities)
    end
end
IndependentBitNoiseModel(probabilitie::Float64, num_bits::Int) = IndependentBitNoiseModel(fill(probabilitie, num_bits))

function to_vector_error_model(em::FactoredBitNoiseModel)
    @assert isvector(em) "The error model must be an vector error model"
    return IndependentBitNoiseModel([em.probabilities[key][2] for key in keys(em.probabilities)])
end

"""`IndependentVectorSampler` is a simple sampler that only works on vector error model."""
struct IndependentVectorSampler <: AbstractSampler end

function sample(em::IndependentBitNoiseModel, num_samples::Int, sampler::IndependentVectorSampler)
    return Array(rand(length(em.probabilities), num_samples) .< em.probabilities)
end

sample(em::FactoredBitNoiseModel, num_samples::Int, sampler::IndependentVectorSampler) = sample(to_vector_error_model(em), num_samples, sampler)

"""
    $TYPEDEF

Represents a quantum error correction decoding problem by combining an error 
model with the code's stabilizer checks and logical operators. This structure 
forms the complete specification needed to perform quantum error correction decoding. 

See also: [`FactoredBitNoiseModel`](@ref)

### Fields
    $TYPEDFIELDS
"""
struct DetectorModelProblem{NMT<:AbstractBitNoiseModel} <: AbstractDecodingProblem
    """Probabilistic error model for the bits"""
    error_model::NMT
    """Parity check matrix
    - **Matrix dimensions**: [`num_checks` × `num_bits`]
    - **Matrix elements**: Boolean (false=0, true=1)
    - **Operation**: Syndrome = `check_matrix` × `error_vector` (mod 2)"""
    check_matrix::Matrix{Bool}
    """Logical operators
    - **Matrix dimensions**: [`num_logicals` × `num_bits`]
    - **Matrix elements**: Boolean (false=0, true=1)"""
    logical_matrix::Matrix{Bool}
    function DetectorModelProblem(error_model::NMT, check_matrix::Matrix{Bool}, logical_matrix::Matrix{Bool}) where NMT<:AbstractBitNoiseModel
        @assert size(check_matrix, 2) == error_model.num_bits "The number of columns of the check matrix must match the number of bits in the error model"
        @assert size(logical_matrix, 2) == error_model.num_bits "The number of columns of the logical matrix must match the number of bits in the error model"
        new{NMT}(error_model, check_matrix, logical_matrix)
    end
end

""" 
    $TYPEDEF

Represents a set of bit string samples. Each column of the following matrices is a sample.

### Fields
    $TYPEDFIELDS
"""
struct BitStringSamples <: AbstractDecodingSamples
    """Physical bits represent the physical state of the bits."""
    physical_bits::Matrix{Bool}
    """Check bits represent the syndrome of the samples."""
    check_bits::Matrix{Bool}
    """Logical bits represent the logical state of samples."""
    logical_bits::Matrix{Bool}
end

function sample(problem::DetectorModelProblem{NMT}, num_samples::Int, sampler::IndependentVectorSampler) where NMT<:FactoredBitNoiseModel
    physical_bits = sample(problem.error_model, num_samples, sampler)
    check_bits = measure_syndrome(problem.check_matrix, physical_bits)
    logical_bits = measure_syndrome(problem.logical_matrix, physical_bits)
    return BitStringSamples(physical_bits, check_bits, logical_bits)
end

function measure_syndrome(check_matrix::AbstractMatrix, error_patterns::AbstractMatrix)
    return Bool.(mod.(check_matrix * error_patterns, 2))
end
measure_syndrome(problem::DetectorModelProblem, error_patterns::Matrix{Bool}) = measure_syndrome(problem.check_matrix, error_patterns)

"""
    $TYPEDEF

Represents a syndrome of a set of bit string samples. Each column of the matrix is a syndrome.

### Fields
    $TYPEDFIELDS
"""
struct MatrixSyndrome <: AbstractSyndrome
    """Check bits represent the syndrome of the samples."""
    check_bits::Matrix{Bool}
end
syndrome(samples::BitStringSamples) = MatrixSyndrome(samples.check_bits)

"""
    $TYPEDEF

Represents a decoding result of a set of bit string samples. Each column of the matrix is a decoded physical state of the bits.

### Fields
    $TYPEDFIELDS
"""
struct MatrixDecodingResult <: AbstractDecodingResult
    """Physical bits represent the decoded physical state of the bits."""
    physical_bits::Matrix{Bool}
end

"""
    $TYPEDSIGNATURES

Calculate the error rate of the decoding result.
"""
function decoding_error_rate(problem::DetectorModelProblem, samples::BitStringSamples, decoding_result::MatrixDecodingResult)
    ep = Bool.(mod.(decoding_result.physical_bits - samples.physical_bits, 2))
    return count(check_decoding_result(problem, ep)) / size(decoding_result.physical_bits, 2)
end

function check_decoding_result(problem::DetectorModelProblem,decoding_result_diff::AbstractMatrix{Bool})
    syndrome_test = any(measure_syndrome(problem.check_matrix, decoding_result_diff), dims=1)
    logical_test = any(measure_syndrome(problem.logical_matrix, decoding_result_diff), dims=1)
    return syndrome_test .|| logical_test
end
