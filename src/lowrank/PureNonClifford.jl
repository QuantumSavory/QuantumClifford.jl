module PureNonClifford

using Random
using Statistics
using LinearAlgebra

import ..QuantumClifford:
    AbstractOperation, AbstractCliffordOperator, AbstractNonCliffordOperator,
    AbstractSingleQubitOperator, AbstractTwoQubitOperator,
    AbstractSymbolicOperator,

    PauliOperator, Stabilizer, MixedDestabilizer,

    apply!, nqubits, stabilizerview, destabilizerview, rank,
    zero, comm,

    sHadamard, sPhase, sX, sCPHASE,

    sMX, sMY, sMZ, sMRX, sMRY, sMRZ, PauliMeasurement,

    SparseGate, CliffordOperator,

    isclifford,

    mctrajectory!, applywstatus!, continue_stat, CircuitStatus,

    sT, sCCZ

using DocStringExtensions

export
    sT, sCCZ,
    LRTrajectoryResults,

    lrtrajectories,
    lrstate,
    lrmeasurements,
    lrcost,

    isclifford,
    stabilizer_extent


"""
$(SIGNATURES)

Return the stabilizer extent ξ(op) for a gate.
For Clifford gates, ξ = 1. Users can extend for custom non-Clifford gates.

The extent determines simulation cost: total cost scales as ∏ⱼ ξ(Vⱼ).

From Proposition 2 of [bravyi2019simulation](@cite): For Clifford magic states, ξ(ψ) = F(ψ)⁻¹
where F(ψ) = max_φ |⟨φ|ψ⟩|² is the stabilizer fidelity.

# Examples
```jldoctest
julia> stabilizer_extent(sHadamard(1))
1.0

julia> round(stabilizer_extent(sT(1)), digits=3)
1.172

julia> stabilizer_extent(sCCZ(1,2,3)) ≈ 16/9
true
```
"""
function stabilizer_extent end

stabilizer_extent(op::AbstractOperation) = isclifford(op) ? 1.0 : error("stabilizer_extent not defined for $(typeof(op)). Please define a method.")

const _T_EXTENT = (cos(π/8) + tan(π/8) * sin(π/8))^2
stabilizer_extent(::sT) = _T_EXTENT

const _CCZ_EXTENT = 16.0 / 9.0
stabilizer_extent(::sCCZ) = _CCZ_EXTENT

"""
$(SIGNATURES)

Sparsification for MixedDestabilizer states (used during incremental simulation).

This is the core routine for incremental sparsification, applying Lemma 6 of [bravyi2019simulation](@cite) to
MixedDestabilizer states that can continue to have gates applied to them.

Samples k ≥ ||c||₁²/δ² states from the decomposition, preserving the L1 norm
of the coefficient vector (in expectation).

# Arguments
- `coefficients::Vector{ComplexF64}`: Current decomposition coefficients
- `states::Vector{MixedDestabilizer}`: Current stabilizer states
- `delta::Float64`: Per-gate approximation error budget
- `rng::AbstractRNG`: Random number generator

# Returns
Named tuple with fields:
- `states::Vector{MixedDestabilizer}`: Sparsified states
- `coefficients::Vector{ComplexF64}`: Sparsified coefficients
- `k::Int`: Number of terms after sparsification
- `l1_norm::Float64`: L1 norm of coefficient vector

# Note
If the current number of terms is already ≤ k, returns copies without sampling.
This avoids unnecessary variance introduction when already sparse enough.
"""
function sparsify_mixed_destabilizer_decomposition(
    coefficients::Vector{ComplexF64}, 
    states::Vector{<:MixedDestabilizer}, 
    delta::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG)
    
    n = length(coefficients)
    n == length(states) || throw(DimensionMismatch(
        "Coefficients ($n) and states ($(length(states))) must have same length"))
    delta > 0 || throw(ArgumentError("δ must be positive, got $delta"))
    n > 0 || throw(ArgumentError("Empty decomposition"))
    
    l1_norm = sum(abs, coefficients)
    l1_norm > 0 || throw(ArgumentError("All coefficients are zero"))
    
    k = max(1, ceil(Int, l1_norm^2 / delta^2))
    
    if n <= k
        return (
            states = [copy(s) for s in states],
            coefficients = copy(coefficients),
            k = n,
            l1_norm = l1_norm
        )
    end
    
    abs_coeffs = abs.(coefficients)
    probabilities = abs_coeffs ./ l1_norm
    cumulative_probs = cumsum(probabilities)
    
    sampled_states = Vector{MixedDestabilizer}(undef, k)
    sampled_coefficients = Vector{ComplexF64}(undef, k)
    uniform_weight = l1_norm / k
    
    for i in 1:k
        r = rand(rng)
        selected_idx = something(findfirst(>=(r), cumulative_probs), n)
        
        phase_factor = coefficients[selected_idx] / abs_coeffs[selected_idx]
        sampled_states[i] = copy(states[selected_idx])
        sampled_coefficients[i] = phase_factor * uniform_weight
    end
    
    return (
        states = sampled_states,
        coefficients = sampled_coefficients,
        k = k,
        l1_norm = l1_norm
    )
end

"""
$(TYPEDEF)

Stabilizer decomposition V|+⟩^⊗t = Σⱼ cⱼ|φⱼ⟩ for Clifford magic states.
Used as intermediate step before applying Lifting Lemma.

For Clifford magic states: ξ(ψ) = F(ψ)⁻¹ where F(ψ) = max_φ |⟨φ|ψ⟩|² 
is the stabilizer fidelity (Proposition 2 of [bravyi2019simulation](@cite)).

$(TYPEDFIELDS)
"""
struct MagicStateDecompositionCache
    "Type of gate (:T, :CCZ, :R_theta)"
    gate_type::Symbol
    "Decomposition coefficients cⱼ"
    coefficients::Vector{ComplexF64}
    "Stabilizer states |φⱼ⟩"
    stabilizer_states::Vector{Stabilizer}
    "||c||₁"
    l1_norm::Float64
    "ξ(ψ) = ||c||₁²"
    stabilizer_extent::Float64
    "F(ψ) = 1/ξ(ψ) for Clifford magic states"
    stabilizer_fidelity::Float64
    
    function MagicStateDecompositionCache(gate_type, coeffs, states)
        length(coeffs) == length(states) || throw(DimensionMismatch(
            "Mismatched coefficients ($(length(coeffs))) and states ($(length(states)))"))
        l1 = sum(abs.(coeffs))
        
        if gate_type == :CCZ
            xi = 16.0/9.0
            fidelity = 9.0/16.0
        elseif gate_type == :R_theta || gate_type == :T
            if length(coeffs) == 2
                c1_real = real(coeffs[1])
                c2_mag = abs(coeffs[2])
                sin_half = c2_mag / sqrt(2)
                cos_half = c1_real + sin_half
                xi = (cos_half + tan(π/8) * sin_half)^2
                fidelity = 1.0 / xi
            else
                xi = l1^2
                fidelity = 1.0 / xi
            end
        else
            xi = l1^2
            fidelity = 1.0 / xi
        end
        
        new(gate_type, coeffs, states, l1, xi, fidelity)
    end
end

"""
$(TYPEDEF)

Sum-over-Cliffords decomposition U = Σⱼ cⱼKⱼ where Kⱼ are Clifford unitaries.
Result of applying Lifting Lemma to magic state decomposition.

$(TYPEDFIELDS)
"""
struct CliffordGateDecompositionCache
    "Type of original non-Clifford gate"
    gate_type::Symbol
    "Decomposition coefficients cⱼ"
    coefficients::Vector{ComplexF64}
    "Clifford circuits Kⱼ"
    clifford_operations::Vector{Vector{AbstractOperation}}
    "||c||₁"
    l1_norm::Float64
    "ξ = ||c||₁²"
    stabilizer_extent::Float64
    "Qubits the gate acts on"
    target_qubits::Vector{Int}
    
    function CliffordGateDecompositionCache(gate_type, coeffs, ops, qubits)
        length(coeffs) == length(ops) || throw(DimensionMismatch(
            "Mismatched coefficients ($(length(coeffs))) and operations ($(length(ops)))"))
        l1 = sum(abs.(coeffs))
        xi = l1^2
        new(gate_type, coeffs, ops, l1, xi, qubits)
    end
end

"""
$(SIGNATURES)

Create magic state decomposition R(θ)|+⟩ = Σⱼ cⱼ|φⱼ⟩ from Eq. (26) of [bravyi2019simulation](@cite).

R(θ)|+⟩ = (cos(θ/2) - sin(θ/2))|+⟩ + √2 sin(θ/2)e^(-iπ/4) S|+⟩

Returns optimal decomposition with ξ(R(θ)) = (cos(θ/2) + tan(π/8)sin(θ/2))².

# Arguments
- `θ::Float64`: Rotation angle
- `nqubits::Int=1`: Number of qubits (for multi-qubit states)

# Returns
`MagicStateDecompositionCache` with the optimal decomposition.
"""
function decompose_rotation_magic_state(θ::Float64; nqubits::Int=1)
    cos_half = cos(θ/2)
    sin_half = sin(θ/2)
    
    c1 = ComplexF64(cos_half - sin_half, 0.0)
    c2 = sqrt(2) * sin_half * exp(-im * π/4)
    
    if nqubits == 1
        plus_state = Stabilizer([PauliOperator(0x00, Bool[true], Bool[false])])
        s_plus_state = Stabilizer([PauliOperator(0x00, Bool[true], Bool[true])])
    else
        plus_generators = [zero(PauliOperator, nqubits) for _ in 1:nqubits]
        for i in 1:nqubits
            plus_generators[i][i] = (true, false)
        end
        plus_state = Stabilizer(plus_generators)
        
        s_plus_generators = [zero(PauliOperator, nqubits) for _ in 1:nqubits]
        for i in 1:nqubits
            s_plus_generators[i][i] = (true, true)
        end
        s_plus_state = Stabilizer(s_plus_generators)
    end
    
    return MagicStateDecompositionCache(:R_theta, [c1, c2], [plus_state, s_plus_state])
end

"""
$(SIGNATURES)

Apply Lifting Lemma (Lemma 1) of [bravyi2019simulation](@cite) to convert R(θ)|+⟩ = Σⱼ cⱼ|φⱼ⟩ to R(θ) = Σⱼ cⱼKⱼ.

For a diagonal gate V, if V|+⟩ = Σⱼ cⱼKⱼ|+⟩ where Kⱼ are diagonal Cliffords,
then V = Σⱼ cⱼKⱼ.

# Arguments
- `magic_decomp`: Magic state decomposition from `decompose_rotation_magic_state`
- `qubit::Int`: Target qubit index

# Returns
`CliffordGateDecompositionCache` with the gate decomposition.
"""
function lifting_lemma_single_qubit(magic_decomp::MagicStateDecompositionCache, qubit::Int)
    coeffs = magic_decomp.coefficients
    
    clifford_ops = Vector{Vector{AbstractOperation}}()
    
    if length(coeffs) == 2
        push!(clifford_ops, AbstractOperation[])
        push!(clifford_ops, [sPhase(qubit)])
    else
        throw(ArgumentError("Unexpected number of terms ($(length(coeffs))) in single-qubit decomposition"))
    end
    
    return CliffordGateDecompositionCache(magic_decomp.gate_type, coeffs, clifford_ops, [qubit])
end

"""
$(SIGNATURES)

Create stabilizer state by applying Clifford operations to |+++⟩.
Internal helper for CCZ decomposition.
"""
function create_ccz_stabilizer_state(operations::Vector{<:AbstractOperation})
    state = MixedDestabilizer(Stabilizer([
        PauliOperator(0x00, Bool[true, false, false], Bool[false, false, false]),
        PauliOperator(0x00, Bool[false, true, false], Bool[false, false, false]),
        PauliOperator(0x00, Bool[false, false, true], Bool[false, false, false])
    ]))
    
    for op in operations
        apply!(state, op)
    end
    
    return Stabilizer(stabilizerview(state))
end

"""
$(SIGNATURES)

Create optimal magic state decomposition for CCZ gate using Proposition 2 of [bravyi2019simulation](@cite).

For Clifford magic states: ξ(ψ) = F(ψ)⁻¹ where F(ψ) = max_φ |⟨φ|ψ⟩|².
For CCZ: F(CCZ) = |⟨+++|CCZ|+++⟩|² = 9/16, so ξ(CCZ) = 16/9.

Uses group decomposition |CCZ⟩ = (1/|Q|⟨CCZ|+++⟩) Σ_{q∈Q} q|+++⟩
where Q = ⟨X₁CZ₂,₃, X₂CZ₁,₃, X₃CZ₁,₂⟩ has 8 elements.

# Returns
`MagicStateDecompositionCache` with optimal decomposition achieving ξ(CCZ) = 16/9.
"""
function decompose_CCZ_magic_state()
    stabilizer_fidelity = 9.0/16.0
    overlap_amplitude = sqrt(stabilizer_fidelity)
    
    group_size = 8
    
    optimal_coeff = 1.0 / (group_size * overlap_amplitude)
    
    states = Vector{Stabilizer}()
    coefficients = ComplexF64[]
    
    group_elements = [
        AbstractOperation[],
        [sX(1), sCPHASE(2,3)],
        [sX(2), sCPHASE(1,3)],
        [sX(3), sCPHASE(1,2)],
        [sX(1), sCPHASE(2,3), sX(2), sCPHASE(1,3)],
        [sX(1), sCPHASE(2,3), sX(3), sCPHASE(1,2)],
        [sX(2), sCPHASE(1,3), sX(3), sCPHASE(1,2)],
        [sX(1), sCPHASE(2,3), sX(2), sCPHASE(1,3), sX(3), sCPHASE(1,2)]
    ]
    
    for operations in group_elements
        state = create_ccz_stabilizer_state(operations)
        push!(states, state)
        push!(coefficients, ComplexF64(optimal_coeff, 0.0))
    end
    
    return MagicStateDecompositionCache(:CCZ, coefficients, states)
end

"""
$(SIGNATURES)

Apply Lifting Lemma to convert CCZ magic state decomposition to gate decomposition.

# Arguments
- `magic_decomp`: Magic state decomposition from `decompose_CCZ_magic_state`
- `qubits::Vector{Int}`: Target qubit indices (must be length 3)

# Returns
`CliffordGateDecompositionCache` with the CCZ gate decomposition.
"""
function lifting_lemma_CCZ(magic_decomp::MagicStateDecompositionCache, qubits::Vector{Int})
    length(qubits) == 3 || throw(ArgumentError("CCZ requires exactly 3 qubits, got $(length(qubits))"))
    
    coeffs = magic_decomp.coefficients
    clifford_ops = Vector{Vector{AbstractOperation}}()
    
    term_operations = [
        AbstractOperation[],
        [sX(qubits[1]), sCPHASE(qubits[2], qubits[3])],
        [sX(qubits[2]), sCPHASE(qubits[1], qubits[3])],
        [sX(qubits[3]), sCPHASE(qubits[1], qubits[2])],
        [sX(qubits[1]), sCPHASE(qubits[2], qubits[3]), sX(qubits[2]), sCPHASE(qubits[1], qubits[3])],
        [sX(qubits[1]), sCPHASE(qubits[2], qubits[3]), sX(qubits[3]), sCPHASE(qubits[1], qubits[2])],
        [sX(qubits[2]), sCPHASE(qubits[1], qubits[3]), sX(qubits[3]), sCPHASE(qubits[1], qubits[2])],
        [sX(qubits[1]), sCPHASE(qubits[2], qubits[3]), sX(qubits[2]), sCPHASE(qubits[1], qubits[3]), sX(qubits[3]), sCPHASE(qubits[1], qubits[2])]
    ]
    
    for ops in term_operations
        push!(clifford_ops, ops)
    end
    
    return CliffordGateDecompositionCache(:CCZ, coeffs, clifford_ops, qubits)
end

"""
$(SIGNATURES)

Get CliffordGateDecompositionCache for a non-Clifford gate.
Uses multiple dispatch - users can extend for custom gates.

# Returns
`CliffordGateDecompositionCache` containing the sum-over-Cliffords decomposition.
"""
function get_gate_decomposition end

function get_gate_decomposition(gate::sT)
    magic_decomp = decompose_rotation_magic_state(π/4; nqubits=1)
    gate_decomp = lifting_lemma_single_qubit(magic_decomp, gate.qubit)
    return CliffordGateDecompositionCache(:T, gate_decomp.coefficients, 
                                          gate_decomp.clifford_operations, [gate.qubit])
end

function get_gate_decomposition(gate::sCCZ)
    magic_decomp = decompose_CCZ_magic_state()
    qubits = collect(gate.qubits)
    return lifting_lemma_CCZ(magic_decomp, qubits)
end

function get_gate_decomposition(gate::AbstractOperation)
    if isclifford(gate)
        throw(ArgumentError("Gate $(typeof(gate)) is Clifford, no decomposition needed"))
    else
        throw(ArgumentError("No decomposition defined for $(typeof(gate)). Please define get_gate_decomposition method."))
    end
end

"""
$(TYPEDEF)

Results from low-rank stabilizer simulation, analogous to PauliFrame results.

$(TYPEDFIELDS)

# Accessing Results
Use `lrmeasurements(result)` to get measurement matrix, similar to `pfmeasurements`.

# Example
```jldoctest
julia> circuit = [sHadamard(1), sT(1), sCNOT(1,2)];

julia> result = lrtrajectories(circuit, 2; trajectories=100, delta=0.1);

julia> measurements = lrmeasurements(result);

julia> size(measurements)
(100, 2)

julia> eltype(measurements)
Bool
```
"""
struct LRTrajectoryResults
    "Measurement outcomes (trajectories × qubits)"
    measurements::Matrix{Bool}
    "Number of sparse stabilizer terms used (k)"
    simulation_cost::Int
    "Target δ parameter"
    approximation_error::Float64
    "Product of gate extents ∏ⱼ ξ(Vⱼ)"
    total_extent::Float64
    "Number of qubits"
    n_qubits::Int
    "Simulation time in seconds"
    total_runtime::Float64
end

function Base.show(io::IO, r::LRTrajectoryResults)
    n_trajectories = size(r.measurements, 1)
    
    println(io, "=== Low-Rank Stabilizer Simulation Results ===")
    println(io, "Qubits: $(r.n_qubits)")
    println(io, "Trajectories: $n_trajectories")
    println(io, "Simulation cost: $(r.simulation_cost) stabilizer terms")
    println(io, "Total extent: $(round(r.total_extent, digits=4))")
    println(io, "Approximation δ: $(r.approximation_error)")
    println(io, "Runtime: $(round(r.total_runtime, digits=2)) seconds")
    
    if n_trajectories > 0
        outcomes = Dict{BitVector, Int}()
        for t in 1:n_trajectories
            bv = BitVector(r.measurements[t, :])
            outcomes[bv] = get(outcomes, bv, 0) + 1
        end
        sorted_outcomes = sort(collect(outcomes), by=x->x[2], rev=true)
        
        println(io, "\nMost frequent outcomes:")
        for (i, (bv, count)) in enumerate(sorted_outcomes[1:min(10, length(sorted_outcomes))])
            bits = join(Int.(bv))
            freq = round(100 * count / n_trajectories, digits=2)
            println(io, "  |$bits⟩: $freq%")
        end
        if length(sorted_outcomes) > 10
            println(io, "  ... and $(length(sorted_outcomes) - 10) other outcomes")
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", r::LRTrajectoryResults)
    show(io, r)
end

"""
$(SIGNATURES)

Extract measurement outcomes from simulation results.
Each row is one trajectory, each column is one measured qubit.

Analogous to `pfmeasurements` for Pauli frames.

# Example
```jldoctest
julia> circuit = [sHadamard(1), sT(1)];

julia> result = lrtrajectories(circuit; trajectories=100);

julia> m = lrmeasurements(result);

julia> size(m)
(100, 1)

julia> 0.0 <= sum(m[:, 1] .== false) / size(m, 1) <= 1.0
true
```
"""
lrmeasurements(r::LRTrajectoryResults) = r.measurements

"""
$(TYPEDEF)
Pure state represented as a weighted sum of stabilizer states:

|ψ⟩ = Σₐ cₐ|φₐ⟩

where:
- `cₐ` are stored in `coefficients`
- `|φₐ⟩` are stored in `states` (as `MixedDestabilizer`)

Similar to `GeneralizedStabilizer` but restricted to pure states. This is the
internal representation used by the low-rank simulation method from [bravyi2019simulation](@cite).

$(TYPEDFIELDS)
"""
mutable struct PureGeneralizedStabilizer
    "MixedDestabilizer states |φₐ⟩ in the decomposition"
    states::Vector{MixedDestabilizer}
    "Coefficients cₐ for each stabilizer state"
    coefficients::Vector{ComplexF64}
    "Per-gate approximation error budget δᵢ"
    delta_per_gate::Float64
    "Running product of gate extents ∏ⱼ ξ(Vⱼ)"
    total_extent::Float64
end

"""
$(SIGNATURES)

Create a `PureGeneralizedStabilizer` initialized to the |0⟩ⁿ state.
"""
function PureGeneralizedStabilizer(n_qubits::Int, delta_per_gate::Float64)
    n_qubits > 0 || throw(ArgumentError("Number of qubits must be positive, got $n_qubits"))
    initial_state = one(MixedDestabilizer, n_qubits)
    PureGeneralizedStabilizer([initial_state], [ComplexF64(1.0)], delta_per_gate, 1.0)
end

nqubits(s::PureGeneralizedStabilizer) = nqubits(s.states[1])

function Base.copy(s::PureGeneralizedStabilizer)
    PureGeneralizedStabilizer(
        [copy(st) for st in s.states],
        copy(s.coefficients),
        s.delta_per_gate,
        s.total_extent
    )
end

"""
$(SIGNATURES)

Apply a Clifford operation to a `PureGeneralizedStabilizer` by applying it to each state.
"""
function apply!(state::PureGeneralizedStabilizer, op::AbstractCliffordOperator)
    for s in state.states
        apply!(s, op)
    end
    state
end

function apply!(state::PureGeneralizedStabilizer, op::AbstractSingleQubitOperator)
    for s in state.states
        apply!(s, op)
    end
    state
end

function apply!(state::PureGeneralizedStabilizer, op::AbstractTwoQubitOperator)
    for s in state.states
        apply!(s, op)
    end
    state
end

function apply!(state::PureGeneralizedStabilizer, op::AbstractSymbolicOperator)
    for s in state.states
        apply!(s, op)
    end
    state
end

function apply!(state::PureGeneralizedStabilizer, op::SparseGate)
    for s in state.states
        apply!(s, op)
    end
    state
end

"""
$(SIGNATURES)

Apply a non-Clifford operation to a `PureGeneralizedStabilizer` by expanding the
sum-over-Cliffords decomposition and then sparsifying.
"""
function apply!(state::PureGeneralizedStabilizer, op::AbstractNonCliffordOperator)
    gate_decomp = get_gate_decomposition(op)
    state.total_extent *= stabilizer_extent(op)

    num_decomp_terms = length(gate_decomp.coefficients)
    current_k = length(state.states)

    expanded_states = Vector{MixedDestabilizer}()
    expanded_coefficients = Vector{ComplexF64}()
    sizehint!(expanded_states, current_k * num_decomp_terms)
    sizehint!(expanded_coefficients, current_k * num_decomp_terms)

    for (s, c) in zip(state.states, state.coefficients)
        for (gc, ops) in zip(gate_decomp.coefficients, gate_decomp.clifford_operations)
            new_s = copy(s)
            for cliff_op in ops
                apply!(new_s, cliff_op)
            end
            push!(expanded_states, new_s)
            push!(expanded_coefficients, c * gc)
        end
    end

    sparse_result = sparsify_mixed_destabilizer_decomposition(
        expanded_coefficients, expanded_states, state.delta_per_gate)
    state.states = sparse_result.states
    state.coefficients = sparse_result.coefficients
    state
end

function applywstatus!(state::PureGeneralizedStabilizer, op::AbstractOperation)
    apply!(state, op), continue_stat
end

"""
$(SIGNATURES)

Create |x⟩ = |x₁x₂...xₙ⟩ as stabilizer state.
Stabilized by (-1)^xᵢ Zᵢ for each qubit i.
"""
function create_computational_basis_state(bitstring::BitVector)
    n = length(bitstring)
    
    generators = [zero(PauliOperator, n) for _ in 1:n]
    for i in 1:n
        generators[i][i] = (false, true)
        if bitstring[i]
            generators[i].phase[] = 0x02
        end
    end
    
    return Stabilizer(generators)
end

"""
$(SIGNATURES)

Compute inner product ⟨state1|state2⟩ using QuantumClifford's dot function.
Based on Section 4.3, Lemma 3 of [bravyi2019simulation](@cite).
"""
function compute_stabilizer_inner_product(state1::Stabilizer, state2::Stabilizer)
    n = nqubits(state1)
    nqubits(state2) == n || throw(DimensionMismatch(
        "States must have same number of qubits: $n vs $(nqubits(state2))"))
    
    return LinearAlgebra.dot(state1, state2)
end

"""
$(SIGNATURES)

Compute amplitude ⟨x|ψ⟩ where |ψ⟩ = Σ cₐ|φₐ⟩ and x is computational basis state.
Uses O(kn³) algorithm from Section 4.3 of [bravyi2019simulation](@cite).
"""
function compute_amplitude(states::Vector{<:Stabilizer}, 
                          coeffs::Vector{ComplexF64},
                          bitstring::BitVector)
    
    basis_state = create_computational_basis_state(bitstring)
    
    amplitude = zero(ComplexF64)
    for (coeff, state) in zip(coeffs, states)
        overlap = compute_stabilizer_inner_product(basis_state, state)
        amplitude += coeff * overlap
    end
    
    return amplitude
end

"""
$(TYPEDEF)

Cached state for efficient Metropolis sampling.
Stores current bitstring and precomputed amplitude to avoid recomputation.

$(TYPEDFIELDS)
"""
mutable struct MetropolisSamplerState
    "Current bitstring being sampled"
    current_bitstring::BitVector
    "Squared amplitude |⟨x|ψ⟩|² for current bitstring"
    current_amplitude_sq::Float64
    "Stabilizer states in the decomposition"
    states::Vector{Stabilizer}
    "Coefficients for each stabilizer state"
    coefficients::Vector{ComplexF64}
    "Number of qubits"
    n_qubits::Int
end

"""
$(SIGNATURES)

Initialize Metropolis sampler with random starting bitstring.
"""
function initialize_metropolis(states::Vector{<:Stabilizer}, 
                              coeffs::Vector{ComplexF64})
    n = nqubits(states[1])
    
    initial_bitstring = BitVector(rand(Bool, n))
    
    amplitude = compute_amplitude(states, coeffs, initial_bitstring)
    amplitude_sq = abs2(amplitude)
    
    return MetropolisSamplerState(initial_bitstring, amplitude_sq, 
                                  Vector{Stabilizer}(states), coeffs, n)
end

"""
$(SIGNATURES)

Perform single Metropolis step: propose bit flip and accept/reject.
Returns true if proposal was accepted.

This is the core O(kn³) operation for each step, where k is the number
of stabilizer states and n³ comes from the inner product computation.
"""
function metropolis_step!(metro::MetropolisSamplerState)
    flip_position = rand(1:metro.n_qubits)
    proposed_bitstring = copy(metro.current_bitstring)
    proposed_bitstring[flip_position] = !proposed_bitstring[flip_position]
    
    proposed_amplitude = compute_amplitude(metro.states, metro.coefficients, proposed_bitstring)
    proposed_amplitude_sq = abs2(proposed_amplitude)
    
    if proposed_amplitude_sq >= metro.current_amplitude_sq
        metro.current_bitstring = proposed_bitstring
        metro.current_amplitude_sq = proposed_amplitude_sq
        return true
    else
        acceptance_ratio = metro.current_amplitude_sq > 0 ? 
            proposed_amplitude_sq / metro.current_amplitude_sq : 1.0
        if rand() < acceptance_ratio
            metro.current_bitstring = proposed_bitstring
            metro.current_amplitude_sq = proposed_amplitude_sq
            return true
        end
    end
    
    return false
end

"""
$(SIGNATURES)

Run burn-in period to reach equilibrium distribution.
Recommended: burn_in ≈ 10n for shallow circuits, 100n for deep circuits.
"""
function metropolis_mixing!(metro::MetropolisSamplerState, burn_in::Int)
    for _ in 1:burn_in
        metropolis_step!(metro)
    end
end

"""
$(SIGNATURES)

MCMC Metropolis sampling from Section 4.2 of [bravyi2019simulation](@cite).
Complexity: O(k × n³ × T) where T = burn_in + n_samples × thinning.

# Arguments
- `state::PureGeneralizedStabilizer`: Sparse stabilizer decomposition from sum-over-Cliffords
- `n_samples::Int`: Number of measurement outcomes to generate
- `burn_in::Int=0`: Mixing time (0 = auto-tune based on circuit)
- `verbose::Bool=false`: Show progress information

# Returns
`Matrix{Bool}` of size (n_samples, n_qubits) with measurement outcomes.
"""
function sample_measurement_outcomes(state::PureGeneralizedStabilizer, 
                                    n_samples::Int; 
                                    burn_in::Int=0,
                                    verbose::Bool=false)
    
    n_samples > 0 || throw(ArgumentError("Number of samples must be positive, got $n_samples"))
    length(state.states) > 0 || throw(ArgumentError("State must contain at least one stabilizer state"))

    stab_states = [Stabilizer(stabilizerview(s)) for s in state.states]
    n_qubits = nqubits(stab_states[1])

    if burn_in == 0
        circuit_depth_proxy = log2(max(1, length(state.states)))
        burn_in = ceil(Int, 10 * n_qubits * (1 + circuit_depth_proxy / 10))
        if verbose
            @info "Auto-tuned burn-in: $burn_in steps"
        end
    end
    
    if verbose
        @info "Initializing Metropolis sampler with $(length(stab_states)) stabilizer states..."
    end
    metro = initialize_metropolis(stab_states, state.coefficients)
    
    if verbose
        @info "Burn-in: $burn_in steps..."
    end
    metropolis_mixing!(metro, burn_in)
    
    measurements = Matrix{Bool}(undef, n_samples, n_qubits)
    acceptance_count = 0
    total_steps = 0
    
    thinning = max(1, div(n_qubits, 2))
    
    if verbose
        @info "Sampling $n_samples outcomes (thinning=$thinning)..."
    end
    
    for s in 1:n_samples
        for _ in 1:thinning
            accepted = metropolis_step!(metro)
            acceptance_count += accepted ? 1 : 0
            total_steps += 1
        end
        
        measurements[s, :] .= metro.current_bitstring
        
        if verbose && s % max(1, div(n_samples, 10)) == 0
            acceptance_rate = acceptance_count / total_steps
            @info "Progress: $s/$n_samples samples (acceptance rate: $(round(acceptance_rate*100, digits=1))%)"
        end
    end
    
    if verbose
        final_acceptance_rate = acceptance_count / total_steps
        @info "Sampling complete: $(round(final_acceptance_rate*100, digits=1))% acceptance rate"
    end
    
    return measurements
end

"""
$(SIGNATURES)

Compute frequency distribution of measurement outcomes.
"""
function compute_outcome_frequencies(measurements::Matrix{Bool})
    n_samples = size(measurements, 1)
    frequency_dict = Dict{BitVector, Float64}()
    
    for s in 1:n_samples
        bv = BitVector(measurements[s, :])
        frequency_dict[bv] = get(frequency_dict, bv, 0.0) + 1.0
    end
    
    for key in keys(frequency_dict)
        frequency_dict[key] /= n_samples
    end
    
    return frequency_dict
end

"""
$(SIGNATURES)

Validation of simulation parameters.
"""
function validate_simulation_parameters(circuit::AbstractVector,
                                       n_qubits::Int,
                                       trajectories::Int,
                                       delta::Float64)
    
    isempty(circuit) && throw(ArgumentError("Circuit cannot be empty"))
    n_qubits > 0 || throw(ArgumentError("Number of qubits must be positive, got $n_qubits"))
    trajectories > 0 || throw(ArgumentError("Number of trajectories must be positive, got $trajectories"))
    0 < delta < 1 || throw(ArgumentError("Delta must be in (0,1), got $delta"))
    
    if trajectories > 100000
        @warn "Large number of trajectories ($trajectories) may result in long simulation times"
    end
    
    if delta < 0.01
        @warn "Very small delta ($delta) may result in large memory usage and slow simulation"
    end
    
    for (idx, op) in enumerate(circuit)
        validate_operation_qubits(op, n_qubits, idx)
    end
    
    return nothing
end

"""
$(SIGNATURES)

Validate that operation uses valid qubit indices.
"""
function validate_operation_qubits(op::AbstractOperation, n_qubits::Int, idx::Int)
    if op isa sT
        1 ≤ op.qubit ≤ n_qubits || throw(ArgumentError(
            "Operation $idx: sT qubit index $(op.qubit) out of range [1,$n_qubits]"))
        return nothing
    elseif op isa sCCZ
        for q in op.qubits
            1 ≤ q ≤ n_qubits || throw(ArgumentError(
                "Operation $idx: sCCZ qubit index $q out of range [1,$n_qubits]"))
        end
        return nothing
    end
    
    op_type = typeof(op)
    
    if hasfield(op_type, :q)
        qubit = op.q
        1 ≤ qubit ≤ n_qubits || throw(ArgumentError(
            "Operation $idx: qubit index $qubit out of range [1,$n_qubits]"))
    end
    
    if hasfield(op_type, :q1) && hasfield(op_type, :q2)
        q1, q2 = op.q1, op.q2
        1 ≤ q1 ≤ n_qubits || throw(ArgumentError(
            "Operation $idx: qubit index $q1 out of range [1,$n_qubits]"))
        1 ≤ q2 ≤ n_qubits || throw(ArgumentError(
            "Operation $idx: qubit index $q2 out of range [1,$n_qubits]"))
        q1 ≠ q2 || throw(ArgumentError(
            "Operation $idx: two-qubit operation cannot target same qubit $q1"))
    end
    
    return nothing
end

"""
$(SIGNATURES)

Infer number of qubits from circuit operations.
"""
function infer_circuit_nqubits(circuit)
    max_qubit = 1
    for op in circuit
        if op isa sT
            max_qubit = max(max_qubit, op.qubit)
        elseif op isa sCCZ
            max_qubit = max(max_qubit, maximum(op.qubits))
        elseif hasfield(typeof(op), :q)
            max_qubit = max(max_qubit, op.q)
        elseif hasfield(typeof(op), :q1) && hasfield(typeof(op), :q2)
            max_qubit = max(max_qubit, op.q1, op.q2)
        end
    end
    return max_qubit
end

"""
$(SIGNATURES)

Simulate quantum circuit with non-Clifford gates using low-rank stabilizer decomposition.

This is the main entry point for non-Clifford simulation, analogous to `pftrajectories` 
for Pauli frame simulation. This method returns **measurement statistics only**, not the 
quantum state. Use [`lrstate`](@ref) if you need access to the state representation.

# Comparison to `mctrajectories` with `GeneralizedStabilizer`
- `lrtrajectories`: Only supports unitary gates (no mid-circuit measurements). 
  Performs implicit Z-basis measurements on all qubits at the end. Faster for 
  circuits with many non-Clifford gates due to sparsification.
- `mctrajectories` with `GeneralizedStabilizer`: Supports mid-circuit measurements 
  and more general operations. Provides access to the full density matrix. More 
  general but slower for large numbers of non-Clifford gates.

# Arguments
- `circuit`: Vector of quantum operations (Clifford gates, sT, sCCZ, etc.)
- `trajectories::Int=1000`: Number of measurement samples to generate
- `delta::Float64=0.1`: Approximation error parameter (smaller = more accurate but slower)
- `verbose::Bool=false`: Show detailed progress information

# Returns
`LRTrajectoryResults` containing measurement outcomes and simulation statistics.
Use `lrmeasurements(result)` to extract the measurement matrix.

# Algorithm
Implements the Sum-over-Cliffords method from Section 2.3.2 of [bravyi2019simulation](@cite)
with incremental sparsification (Section 5.2) after each non-Clifford gate:

1. Initialize with |0⟩ state
2. For each gate:
   - Clifford: apply to all current stabilizer states
   - Non-Clifford: expand decomposition, then sparsify immediately
3. Sample measurement outcomes via Metropolis MCMC (Section 4.2) with implicit Z-basis measurement

# Performance
- Memory: O(k × n) where k ≈ total_extent/δ² (never stores 2^m terms)
- Simulation cost scales as ∏ⱼ ξ(Vⱼ) where Vⱼ are the non-Clifford gates and ξ(T) ≈ 1.172, ξ(CCZ) = 16/9
- Time complexity: O(m × k² × n³) for simulation + O(k × n³ × trajectories) for sampling

# Example
```jldoctest
julia> circuit = [sHadamard(1), sT(1), sHadamard(1)];

julia> result = lrtrajectories(circuit; trajectories=100, delta=0.1);

julia> result.simulation_cost > 0
true

julia> measurements = lrmeasurements(result);

julia> size(measurements)
(100, 1)

julia> p0 = sum(measurements[:, 1] .== false) / size(measurements, 1);

julia> 0.0 <= p0 <= 1.0
true
```

See also: [`lrstate`](@ref), [`lrmeasurements`](@ref), [`lrcost`](@ref)
"""
function lrtrajectories(circuit::AbstractVector{<:AbstractOperation};
                        trajectories::Int=1000,
                        delta::Float64=0.1,
                        verbose::Bool=false)

    n_qubits = infer_circuit_nqubits(circuit)
    return lrtrajectories(circuit, n_qubits; trajectories, delta, verbose)
end

function lrtrajectories(circuit::AbstractVector{<:AbstractOperation},
                        n_qubits::Int;
                        trajectories::Int=1000,
                        delta::Float64=0.1,
                        verbose::Bool=false)

    validate_simulation_parameters(circuit, n_qubits, trajectories, delta)

    start_time = time()

    state = _lrstate(circuit, n_qubits, delta)
    measurements = sample_measurement_outcomes(state, trajectories; verbose)
    total_runtime = time() - start_time

    return LRTrajectoryResults(
        measurements,
        length(state.states),
        delta,
        state.total_extent,
        n_qubits,
        total_runtime
    )
end

function lrtrajectories(circuit::AbstractVector, n_qubits::Int; kwargs...)
    isempty(circuit) && throw(ArgumentError("Circuit cannot be empty"))

    for (i, op) in enumerate(circuit)
        if !(op isa AbstractOperation)
            throw(ArgumentError("Element $i is not an AbstractOperation: got $(typeof(op))"))
        end
    end

    typed_circuit = AbstractOperation[op for op in circuit]
    return lrtrajectories(typed_circuit, n_qubits; kwargs...)
end

"""
$(SIGNATURES)

Run low-rank simulation and return the state representation.

Unlike [`lrtrajectories`](@ref) which returns measurement samples, this function returns
the `PureGeneralizedStabilizer` state directly. This is useful when you need
the state for further analysis or custom sampling.

# Arguments
- `circuit`: Vector of quantum operations (Clifford gates, sT, sCCZ, etc.)
- `n_qubits::Int`: Number of qubits (optional, inferred from circuit if not provided)
- `delta::Float64=0.1`: Approximation error parameter
- `verbose::Bool=false`: Show detailed progress information

# Returns
`PureGeneralizedStabilizer` containing the sparse stabilizer decomposition |ψ⟩ = Σₐ cₐ|φₐ⟩.

# Example
```jldoctest
julia> circuit = [sHadamard(1), sT(1)];

julia> state = lrstate(circuit; delta=0.1);

julia> length(state.coefficients) > 0
true
```

See also: [`lrtrajectories`](@ref), [`PureGeneralizedStabilizer`](@ref)
"""
function lrstate(circuit::AbstractVector{<:AbstractOperation};
                 delta::Float64=0.1,
                 verbose::Bool=false)
    n_qubits = infer_circuit_nqubits(circuit)
    return lrstate(circuit, n_qubits; delta, verbose)
end

function lrstate(circuit::AbstractVector{<:AbstractOperation},
                 n_qubits::Int;
                 delta::Float64=0.1,
                 verbose::Bool=false)
    return _lrstate(circuit, n_qubits, delta)
end

"""
$(SIGNATURES)

Internal helper: create PureGeneralizedStabilizer and run circuit via mctrajectory!.
"""
function _lrstate(circuit::AbstractVector{<:AbstractOperation}, n_qubits::Int, delta::Float64)
    non_clifford_count = count(op -> !isclifford(op), circuit)
    delta_per_gate = non_clifford_count > 0 ? delta / non_clifford_count : delta

    state = PureGeneralizedStabilizer(n_qubits, delta_per_gate)
    mctrajectory!(state, circuit)
    state
end

"""
$(SIGNATURES)

Estimate simulation cost before running full simulation.

This is useful for determining if a circuit is feasible to simulate
and for tuning the delta parameter.

# Arguments
- `circuit`: Vector of quantum operations
- `delta::Float64=0.1`: Approximation error parameter

# Returns
Named tuple with:
- `estimated_k::Int`: Number of stabilizer terms after sparsification (k ≤ 1 + ξ/δ²)
- `total_extent::Float64`: Product ∏ⱼ ξ(Vⱼ) of gate extents
- `non_clifford_count::Int`: Number of non-Clifford gates
- `gate_extents::Vector{Float64}`: Individual gate extents
- `scaling_factor::Float64`: Geometric mean of gate extents

# Performance Guidelines
- k < 1000: Fast simulation (seconds)
- k < 100,000: Moderate simulation (minutes)
- k > 1,000,000: Long simulation (hours), consider larger δ

# Note on Incremental Sparsification
With incremental sparsification, actual memory usage during simulation
is bounded by the maximum k at any single gate, not the final k.
This is typically much smaller than the estimated final k.

# Example
```jldoctest
julia> circuit = [sHadamard(1), sT(1), sT(1), sT(1), sHadamard(1)];

julia> cost = lrcost(circuit; delta=0.1);

julia> cost.non_clifford_count
3

julia> cost.estimated_k > 0
true

julia> cost.total_extent ≈ stabilizer_extent(sT(1))^3
true
```
"""
function lrcost(circuit::AbstractVector{<:AbstractOperation}; delta::Float64=0.1)
    total_extent = 1.0
    non_clifford_count = 0
    gate_extents = Float64[]
    
    for op in circuit
        if !isclifford(op)
            xi = stabilizer_extent(op)
            total_extent *= xi
            non_clifford_count += 1
            push!(gate_extents, xi)
        end
    end
    
    estimated_k = if non_clifford_count == 0
        1
    else
        ceil(Int, 1 + total_extent / delta^2)
    end
    
    scaling_factor = non_clifford_count > 0 ? 
        total_extent^(1/non_clifford_count) : 1.0
    
    return (
        estimated_k = estimated_k,
        total_extent = total_extent,
        non_clifford_count = non_clifford_count,
        gate_extents = gate_extents,
        scaling_factor = scaling_factor
    )
end

end