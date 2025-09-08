"""
    $TYPEDEF

Represents a probabilistic error model for a system of bits using a factored 
distribution approach. The model decomposes the joint probability distribution 
over all bits into smaller, manageable distributions over subsets of bits.

### Fields
    $TYPEDFIELDS
"""
struct ErrorModel{VT<:AbstractArray{Float64}}
    """The number of bits in the error model."""
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
    function ErrorModel(num_bits::Int, probabilities::Dict{Vector{Int},VT}) where VT<:AbstractArray{Float64}
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
function ErrorModel(probabilities::Vector{Float64})
    num_bits = length(probabilities)
    probabilities = Dict([i] => [1.0 - probabilities[i], probabilities[i]] for i in 1:num_bits)
    return ErrorModel(num_bits, probabilities)
end

"""
    $TYPEDSIGNATURES

Create a depolarization error model from a vector of error probabilities. 
The `i`-th element of the vector is the depolarization probability of qubit `i`.
"""
function depolarization_error_model(pvec::Vector{Float64}, qubit_num::Int)
    @assert length(pvec) == qubit_num "The length of the vector of error probabilities must match the number of qubits"
    probabilities = Dict([i, i + qubit_num] => [1.0-pvec[i] pvec[i]/3; pvec[i]/3 pvec[i]/3] for i in 1:qubit_num)
    return ErrorModel(2 * qubit_num, probabilities)
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
isvector(em::ErrorModel) = all(length(key) == 1 for key in keys(em.probabilities))

"""
    $TYPEDSIGNATURES

Check if the error model is an independent error model. Independent error model means that 
"""
function isindependent(em::ErrorModel)
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

abstract type AbstractSampler end
"""`IndependentVectorSampler` is a simple sampler that only works on vector error model."""
struct IndependentVectorSampler <: AbstractSampler end

"""
    $TYPEDSIGNATURES

Generate a random error pattern from the error model. `sampler` is the sampler to use.

See also: [`IndependentVectorSampler`](@ref)
"""
function random_error_pattern(em::ErrorModel, num_samples::Int, sampler::IndependentVectorSampler)
    @assert isvector(em) "The error model must be a vector error model"
    return [rand() < em.probabilities[[i]][2] for i in 1:em.num_bits, _ in 1:num_samples]
end

"""
    $TYPEDEF

Represents a quantum error correction decoding problem by combining an error 
model with the code's stabilizer checks and logical operators. This structure 
forms the complete specification needed to perform quantum error correction decoding. 

See also: [`ErrorModel`](@ref)

### Fields
    $TYPEDFIELDS
"""
struct DecodingProblem{VT}
    """Probabilistic error model for the bits"""
    error_model::ErrorModel{VT}
    """Parity check matrix
    - **Matrix dimensions**: [num_checks × num_bits]
    - **Matrix elements**: Boolean (false=0, true=1)
    - **Operation**: Syndrome = check_matrix × error_vector (mod 2)"""
    check_matrix::Matrix{Bool}
    """Logical operators
    - **Matrix dimensions**: [num_logicals × num_bits]
    - **Matrix elements**: Boolean (false=0, true=1)"""
    logical_matrix::Matrix{Bool}
    function DecodingProblem(error_model::ErrorModel{VT}, check_matrix::Matrix{Bool}, logical_matrix::Matrix{Bool}) where VT
        @assert size(check_matrix, 2) == error_model.num_bits "The number of columns of the check matrix must match the number of bits in the error model"
        @assert size(logical_matrix, 2) == error_model.num_bits "The number of columns of the logical matrix must match the number of bits in the error model"
        new{VT}(error_model, check_matrix, logical_matrix)
    end
end
# TODO: generate code capacity noise decoding problem. Need to compute the logical operators first.

"""
    $TYPEDSIGNATURES

Extract the syndrome from the error patterns.
"""
function syndrome_extraction(check_matrix::AbstractMatrix, error_patterns::AbstractMatrix)
    return Bool.(mod.(check_matrix * error_patterns, 2))
end
syndrome_extraction(problem::DecodingProblem, error_patterns::Matrix{Bool}) = syndrome_extraction(problem.check_matrix, error_patterns)


abstract type AbstractDecoder end

"""
    decode(problem::DecodingProblem,syndrome::Matrix{Bool},decoder::AbstractDecoder)

Decode the error pattern from the syndrome.
"""
function decode end

"""
    $TYPEDSIGNATURES

Check if the error pattern is a logical error.
"""
function check_decoding_result(error_pattern1::AbstractMatrix, error_pattern2::AbstractMatrix, problem::DecodingProblem)
    ep = error_pattern1 - error_pattern2
    syndrome_test = any(syndrome_extraction(problem.check_matrix, ep), dims=1)
    logical_test = any(syndrome_extraction(problem.logical_matrix, ep), dims=1)
    return syndrome_test .|| logical_test
end
