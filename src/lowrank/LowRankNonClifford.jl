module LowRankNonClifford

using Random
using Statistics
using LinearAlgebra

import ..QuantumClifford:
    AbstractOperation, AbstractCliffordOperator,
    AbstractSingleQubitOperator, AbstractTwoQubitOperator,
    AbstractSymbolicOperator,
    
    PauliOperator, Stabilizer, MixedDestabilizer,
    
    apply!, nqubits, stabilizerview, destabilizerview, rank,
    zero, comm,
    
    sHadamard, sPhase, sX, sCPHASE,
    
    sMX, sMY, sMZ, sMRX, sMRY, sMRZ, PauliMeasurement,
    
    SparseGate, CliffordOperator

using DocStringExtensions

export
    TGate, CCZGate,
    LRTrajectoryResults,
    
    lrtrajectories,
    lrmeasurements,
    lrcost,
    
    isclifford,
    stabilizer_extent

"""
    isclifford(op::AbstractOperation) -> Bool

Trait function to determine if an operation is a Clifford gate.
Users can extend this for custom gate types by defining new methods.

This replaces the previous `identify_gate_type` function with an extensible
multimethod approach following Julia idioms.

# Examples
```julia
isclifford(sHadamard(1))   # true
isclifford(sPhase(1))      # true
isclifford(sCNOT(1,2))     # true
isclifford(TGate(1))       # false
isclifford(CCZGate(1,2,3)) # false

# Extend for custom gates:
isclifford(::MyCustomCliffordGate) = true
isclifford(::MyNonCliffordGate) = false
```

"""
function isclifford end

isclifford(::AbstractCliffordOperator) = true
isclifford(::AbstractSingleQubitOperator) = true
isclifford(::AbstractTwoQubitOperator) = true
isclifford(::AbstractSymbolicOperator) = true

isclifford(::sMX) = true
isclifford(::sMY) = true
isclifford(::sMZ) = true
isclifford(::sMRX) = true
isclifford(::sMRY) = true
isclifford(::sMRZ) = true
isclifford(::PauliMeasurement) = true

isclifford(op::SparseGate) = isclifford(op.cliff)

isclifford(::AbstractOperation) = false

"""
    stabilizer_extent(op::AbstractOperation) -> Float64

Return the stabilizer extent ξ(op) for a gate.
For Clifford gates, ξ = 1. Users can extend for custom non-Clifford gates.

The extent determines simulation cost: total cost scales as ∏ⱼ ξ(Vⱼ).

From Proposition 2 of Bravyi et al.: For Clifford magic states, ξ(ψ) = F(ψ)⁻¹
where F(ψ) = max_φ |⟨φ|ψ⟩|² is the stabilizer fidelity.

# Examples
```julia
stabilizer_extent(sHadamard(1))   # 1.0 (Clifford)
stabilizer_extent(TGate(1))       # ≈ 1.172
stabilizer_extent(CCZGate(1,2,3)) # 16/9 ≈ 1.778

# Extend for custom gates:
stabilizer_extent(::MyGate) = 2.5
```

"""
function stabilizer_extent end

stabilizer_extent(op::AbstractOperation) = isclifford(op) ? 1.0 : error("stabilizer_extent not defined for $(typeof(op)). Please define a method.")

"""
    TGate <: AbstractOperation

T gate (π/8 phase rotation). This is a non-Clifford gate with optimal 
stabilizer extent ξ(T) = (cos(π/8) + tan(π/8)sin(π/8))² ≈ 1.172.

The T gate is diagonal in the computational basis:
T = diag(1, e^(iπ/4)) = R_z(π/4)

# Fields
- `qubit::Int`: Target qubit index (1-indexed)

# Example
```julia
circuit = [sHadamard(1), TGate(1), sHadamard(1)]
result = lrtrajectories(circuit, 1; trajectories=1000)
```

"""
struct TGate <: AbstractOperation
    qubit::Int
    
    function TGate(q::Int)
        q > 0 || throw(ArgumentError("Qubit index must be positive, got $q"))
        new(q)
    end
end

nqubits(::TGate) = 1
isclifford(::TGate) = false

const _T_EXTENT = (cos(π/8) + tan(π/8) * sin(π/8))^2
stabilizer_extent(::TGate) = _T_EXTENT

"""
    CCZGate <: AbstractOperation

Controlled-Controlled-Z gate. This is a non-Clifford gate with optimal 
stabilizer extent ξ(CCZ) = 16/9 ≈ 1.778.

The CCZ gate applies a phase of -1 when all three qubits are in state |1⟩:
CCZ|xyz⟩ = (-1)^(xyz)|xyz⟩

# Fields
- `qubits::NTuple{3,Int}`: Target qubit indices (sorted, 1-indexed)

# Example
```julia
circuit = [sHadamard(1), sHadamard(2), sHadamard(3), CCZGate(1, 2, 3)]
result = lrtrajectories(circuit, 3; trajectories=1000)
```

"""
struct CCZGate <: AbstractOperation
    qubits::NTuple{3, Int}
    
    function CCZGate(q1::Int, q2::Int, q3::Int)
        all(q -> q > 0, (q1, q2, q3)) || throw(ArgumentError("All qubit indices must be positive"))
        length(unique((q1, q2, q3))) == 3 || throw(ArgumentError("CCZ requires 3 distinct qubits, got $q1, $q2, $q3"))
        new(Tuple(sort([q1, q2, q3])))
    end
end

CCZGate(qubits::Vector{Int}) = CCZGate(qubits[1], qubits[2], qubits[3])

nqubits(::CCZGate) = 3
isclifford(::CCZGate) = false

const _CCZ_EXTENT = 16.0 / 9.0
stabilizer_extent(::CCZGate) = _CCZ_EXTENT

"""
    SparsifiedDecomposition

Result of applying Sparsification Lemma (Lemma 6) from Section 5.2.
Represents |Ω⟩ = (||c||₁/k) Σₐ₌₁ᵏ |ωₐ⟩ that approximates |ψ⟩ = Σⱼ cⱼ|φⱼ⟩.

# Fields
- `states::Vector{Stabilizer}`: Sampled stabilizer states |ωₐ⟩
- `coefficients::Vector{ComplexF64}`: Corresponding coefficients
- `k::Int`: Number of terms in sparse approximation
- `original_l1_norm::Float64`: ||c||₁ of original decomposition
- `approximation_error_bound::Float64`: Target δ parameter
- `expected_norm_bound::Float64`: E[||Ω||²] = 1 + ||c||₁²/k

"""
struct SparsifiedDecomposition
    states::Vector{Stabilizer}
    coefficients::Vector{ComplexF64}
    k::Int
    original_l1_norm::Float64
    approximation_error_bound::Float64
    expected_norm_bound::Float64
    
    function SparsifiedDecomposition(states, coeffs, k, l1_norm, delta)
        length(states) == length(coeffs) == k || throw(DimensionMismatch(
            "Inconsistent array lengths: states=$(length(states)), coeffs=$(length(coeffs)), k=$k"))
        expected_norm = 1.0 + l1_norm^2 / k
        new(states, coeffs, k, l1_norm, delta, expected_norm)
    end
end

"""
    sparsify_stabilizer_decomposition(coefficients, states, delta; rng=Random.GLOBAL_RNG)

Implements Sparsification Lemma (Lemma 6) from Section 5.2 of Bravyi et al.

Given |ψ⟩ = Σⱼ cⱼ|φⱼ⟩ with L1 norm ||c||₁, constructs random state 
|Ω⟩ = (||c||₁/k) Σₐ₌₁ᵏ |ωₐ⟩ where each |ωₐ⟩ is sampled from {|φⱼ⟩} 
with probability |cⱼ|/||c||₁.

The approximation satisfies E[||ψ - Ω||²] = ||c||₁²/k ≤ δ² for k ≥ ||c||₁²/δ².

# Arguments
- `coefficients::Vector{ComplexF64}`: Coefficients cⱼ in decomposition
- `states::Vector{<:Stabilizer}`: Stabilizer states |φⱼ⟩
- `delta::Float64`: Target approximation error (must be positive)
- `rng::AbstractRNG`: Random number generator (default: global RNG)

# Returns
`SparsifiedDecomposition` containing the sparse approximation.

# Mathematical Guarantee (Theorem 1)
For the returned sparse state: χ_δ(ψ) ≤ 1 + ||c||₁²/δ²

"""
function sparsify_stabilizer_decomposition(coefficients::Vector{ComplexF64}, 
                                          states::Vector{<:Stabilizer}, 
                                          delta::Float64;
                                          rng::AbstractRNG=Random.GLOBAL_RNG)
    
    length(coefficients) == length(states) || throw(DimensionMismatch(
        "Coefficients and states must have same length: $(length(coefficients)) vs $(length(states))"))
    delta > 0 || throw(ArgumentError("Approximation error δ must be positive, got $delta"))
    
    l1_norm = sum(abs.(coefficients))
    l1_norm > 0 || throw(ArgumentError("All coefficients are zero - no valid decomposition"))
    
    k = max(1, ceil(Int, l1_norm^2 / delta^2))
    
    abs_coeffs = abs.(coefficients)
    probabilities = abs_coeffs / l1_norm
    cumulative_probs = cumsum(probabilities)
    
    sampled_states = Vector{Stabilizer}()
    sampled_coefficients = ComplexF64[]
    
    uniform_weight = l1_norm / k
    
    for i in 1:k
        r = rand(rng)
        selected_idx = findfirst(p -> p >= r, cumulative_probs)
        if selected_idx === nothing
            selected_idx = length(states)
        end
        
        phase_factor = coefficients[selected_idx] / abs_coeffs[selected_idx]
        push!(sampled_states, copy(states[selected_idx]))
        push!(sampled_coefficients, phase_factor * uniform_weight)
    end
    
    return SparsifiedDecomposition(sampled_states, sampled_coefficients, k, l1_norm, delta)
end

"""
    sparsify_mixed_destabilizer_decomposition(coefficients, states, delta; rng)

Sparsification for MixedDestabilizer states (used during incremental simulation).

This is the core routine for incremental sparsification, applying Lemma 6 to
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
    estimate_sparsification_quality(sparse::SparsifiedDecomposition)

Estimate quality bounds from Lemma 7 (Sparsification tail bound).

# Returns
Named tuple with:
- `k`: Number of sparse terms
- `expected_error`: E[||ψ - Ω||²] = ||c||₁²/k
- `error_bound`: Target δ parameter
- `expected_norm`: E[||Ω||²] = 1 + ||c||₁²/k

"""
function estimate_sparsification_quality(sparse::SparsifiedDecomposition)
    return (
        k = sparse.k,
        expected_error = sparse.original_l1_norm^2 / sparse.k,
        error_bound = sparse.approximation_error_bound,
        expected_norm = sparse.expected_norm_bound
    )
end

"""
    MagicStateDecompositionCache

Stabilizer decomposition V|+⟩^⊗t = Σⱼ cⱼ|φⱼ⟩ for Clifford magic states.
Used as intermediate step before applying Lifting Lemma.

For Clifford magic states: ξ(ψ) = F(ψ)⁻¹ where F(ψ) = max_φ |⟨φ|ψ⟩|² 
is the stabilizer fidelity (Proposition 2).

# Fields
- `gate_type::Symbol`: Type of gate (:T, :CCZ, :R_theta)
- `coefficients::Vector{ComplexF64}`: Decomposition coefficients cⱼ
- `stabilizer_states::Vector{Stabilizer}`: Stabilizer states |φⱼ⟩
- `l1_norm::Float64`: ||c||₁
- `stabilizer_extent::Float64`: ξ(ψ) = ||c||₁²
- `stabilizer_fidelity::Float64`: F(ψ) = 1/ξ(ψ) for Clifford magic states

"""
struct MagicStateDecompositionCache
    gate_type::Symbol
    coefficients::Vector{ComplexF64}
    stabilizer_states::Vector{Stabilizer}
    l1_norm::Float64
    stabilizer_extent::Float64
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
    CliffordGateDecompositionCache

Sum-over-Cliffords decomposition U = Σⱼ cⱼKⱼ where Kⱼ are Clifford unitaries.
Result of applying Lifting Lemma to magic state decomposition.

# Fields
- `gate_type::Symbol`: Type of original non-Clifford gate
- `coefficients::Vector{ComplexF64}`: Decomposition coefficients cⱼ
- `clifford_operations::Vector{Vector{AbstractOperation}}`: Clifford circuits Kⱼ
- `l1_norm::Float64`: ||c||₁
- `stabilizer_extent::Float64`: ξ = ||c||₁²
- `target_qubits::Vector{Int}`: Qubits the gate acts on
"""
struct CliffordGateDecompositionCache
    gate_type::Symbol
    coefficients::Vector{ComplexF64}
    clifford_operations::Vector{Vector{AbstractOperation}}
    l1_norm::Float64
    stabilizer_extent::Float64
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
    decompose_rotation_magic_state(θ; nqubits=1)

Create magic state decomposition R(θ)|+⟩ = Σⱼ cⱼ|φⱼ⟩ from Eq. (26).

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
    lifting_lemma_single_qubit(magic_decomp::MagicStateDecompositionCache, qubit::Int)

Apply Lifting Lemma (Lemma 1) to convert R(θ)|+⟩ = Σⱼ cⱼ|φⱼ⟩ to R(θ) = Σⱼ cⱼKⱼ.

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
    decompose_T_gate(qubit::Int)

T gate decomposition: T = R(π/4) with ξ(T) = (cos(π/8) + tan(π/8)sin(π/8))² ≈ 1.172.

# Arguments
- `qubit::Int`: Target qubit index

# Returns
`CliffordGateDecompositionCache` with optimal T gate decomposition.

"""
function decompose_T_gate(qubit::Int)
    magic_decomp = decompose_rotation_magic_state(π/4; nqubits=1)
    gate_decomp = lifting_lemma_single_qubit(magic_decomp, qubit)
    return CliffordGateDecompositionCache(:T, gate_decomp.coefficients, 
                                          gate_decomp.clifford_operations, [qubit])
end

"""
    create_ccz_stabilizer_state(operations::Vector{<:AbstractOperation})

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
    decompose_CCZ_magic_state()

Create optimal magic state decomposition for CCZ gate using Proposition 2.

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
    lifting_lemma_CCZ(magic_decomp::MagicStateDecompositionCache, qubits::Vector{Int})

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
    decompose_CCZ_gate(qubits::Vector{Int})

Get optimal CCZ gate decomposition with ξ(CCZ) = 16/9.

# Arguments
- `qubits::Vector{Int}`: Target qubit indices (must be length 3)

# Returns
`CliffordGateDecompositionCache` with optimal CCZ gate decomposition.

"""
function decompose_CCZ_gate(qubits::Vector{Int})
    magic_decomp = decompose_CCZ_magic_state()
    return lifting_lemma_CCZ(magic_decomp, qubits)
end

decompose_CCZ_gate(q1::Int, q2::Int, q3::Int) = decompose_CCZ_gate([q1, q2, q3])

"""
    get_gate_decomposition(gate::AbstractOperation)

Get CliffordGateDecompositionCache for a non-Clifford gate.
Uses multiple dispatch - users can extend for custom gates.

# Returns
`CliffordGateDecompositionCache` containing the sum-over-Cliffords decomposition.
"""
function get_gate_decomposition end

function get_gate_decomposition(gate::TGate)
    return decompose_T_gate(gate.qubit)
end

function get_gate_decomposition(gate::CCZGate)
    return decompose_CCZ_gate(collect(gate.qubits))
end

function get_gate_decomposition(gate::AbstractOperation)
    if isclifford(gate)
        throw(ArgumentError("Gate $(typeof(gate)) is Clifford, no decomposition needed"))
    else
        throw(ArgumentError("No decomposition defined for $(typeof(gate)). Please define get_gate_decomposition method."))
    end
end

"""
    LRTrajectoryResults

Results from low-rank stabilizer simulation, analogous to PauliFrame results.

# Fields
- `measurements::Matrix{Bool}`: Measurement outcomes (trajectories × qubits)
- `simulation_cost::Int`: Number of sparse stabilizer terms used (k)
- `approximation_error::Float64`: Target δ parameter
- `total_extent::Float64`: Product of gate extents ∏ⱼ ξ(Vⱼ)
- `n_qubits::Int`: Number of qubits
- `total_runtime::Float64`: Simulation time in seconds

# Accessing Results
Use `lrmeasurements(result)` to get measurement matrix, similar to `pfmeasurements`.

# Example
```julia
result = lrtrajectories(circuit, 2; trajectories=1000, delta=0.1)
measurements = lrmeasurements(result)  # Matrix{Bool} of size (1000, 2)
```

"""
struct LRTrajectoryResults
    measurements::Matrix{Bool}
    simulation_cost::Int
    approximation_error::Float64
    total_extent::Float64
    n_qubits::Int
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
    lrmeasurements(result::LRTrajectoryResults) -> Matrix{Bool}

Extract measurement outcomes from simulation results.
Each row is one trajectory, each column is one measured qubit.

Analogous to `pfmeasurements` for Pauli frames.

# Example
```julia
result = lrtrajectories(circuit; trajectories=1000)
m = lrmeasurements(result)
prob_0 = sum(m[:, 1] .== false) / size(m, 1)  # Probability of |0⟩ on qubit 1
```

"""
lrmeasurements(r::LRTrajectoryResults) = r.measurements

"""
    SimulationState

Internal state after sum-over-Cliffords simulation completes.
Contains the sparse stabilizer decomposition ready for sampling.
"""
struct SimulationState
    sparse_states::Vector{Stabilizer}
    coefficients::Vector{ComplexF64}
    simulation_cost::Int
    approximation_error::Float64
    original_extent::Float64
end

"""
    create_computational_zero_state(n_qubits::Int) -> MixedDestabilizer

Creates proper |0ⁿ⟩ state stabilized by Z₁, Z₂, ..., Zₙ.
Used as initial state for Sum-over-Cliffords method.
"""
function create_computational_zero_state(n_qubits::Int)
    n_qubits > 0 || throw(ArgumentError("Number of qubits must be positive, got $n_qubits"))
    
    z_operators = [zero(PauliOperator, n_qubits) for _ in 1:n_qubits]
    for i in 1:n_qubits
        z_operators[i][i] = (false, true)
    end
    
    return MixedDestabilizer(Stabilizer(z_operators))
end

"""
    create_computational_basis_state(bitstring::BitVector) -> Stabilizer

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
    compute_stabilizer_inner_product(state1::Stabilizer, state2::Stabilizer) -> ComplexF64

Compute inner product ⟨state1|state2⟩ using QuantumClifford's dot function.
Based on Section 4.3, Lemma 3 of Bravyi et al. 2019.
"""
function compute_stabilizer_inner_product(state1::Stabilizer, state2::Stabilizer)
    n = nqubits(state1)
    nqubits(state2) == n || throw(DimensionMismatch(
        "States must have same number of qubits: $n vs $(nqubits(state2))"))
    
    return LinearAlgebra.dot(state1, state2)
end

"""
    compute_amplitude(states, coeffs, bitstring) -> ComplexF64

Compute amplitude ⟨x|ψ⟩ where |ψ⟩ = Σ cₐ|φₐ⟩ and x is computational basis state.
Uses O(kn³) algorithm from Section 4.3.
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
    simulate_sum_over_cliffords(circuit, n_qubits, delta; verbose=false, rng=GLOBAL_RNG)

Incremental Sum-over-Cliffords simulation with per-gate sparsification.

This prevents exponential memory blowup
by applying the Sparsification Lemma (Lemma 6) after each non-Clifford gate,
rather than building the full 2^m decomposition first

# Algorithm (Section 2.3.2 + Section 5.2)
For each gate in circuit:
1. If Clifford: apply to all current stabilizer states in-place
2. If non-Clifford:
   a. Get gate decomposition U = Σⱼ cⱼKⱼ (e.g., T = c₁I + c₂S)
   b. Expand: each current state |φₐ⟩ becomes {Kⱼ|φₐ⟩} with coefficients multiplied
   c. Sparsify immediately using Lemma 6 with per-gate error budget δᵢ

# Error Budget Management
Total approximation error δ is distributed as δᵢ = δ/m per gate, where m is
the number of non-Clifford gates. By triangle inequality on L2 norms:
  ||ψ_exact - ψ_approx|| ≤ Σᵢ δᵢ = δ

# Complexity Analysis
- Memory: O(k_max × n) where k_max = max over gates of (||c||₁²/δᵢ²)
  For T-gates: k_max ≈ ξ(T)^(gates_so_far) / δᵢ² 
- Time: O(m × k² × decomp_size × n³) for gate application and inner products
- Contrast with non-incremental: would require O(2^m × n) memory

# Mathematical Guarantees
- L1 norm preserved through sparsification (by construction of Lemma 6)
- After m gates: ||c||₁ = √(ξ₁ × ξ₂ × ... × ξₘ) = √(total_extent)
- Final k ≤ 1 + total_extent/δ² (matches Theorem 1)

# Arguments
- `circuit::Vector{<:AbstractOperation}`: Quantum circuit to simulate
- `n_qubits::Int`: Number of qubits
- `delta::Float64`: Total approximation error bound
- `verbose::Bool=false`: Print progress information
- `rng::AbstractRNG=Random.GLOBAL_RNG`: Random number generator for sparsification

# Returns
`SimulationState` containing sparse stabilizer decomposition ready for sampling.

# References
- Section 2.3.2: Sum-over-Cliffords method
- Section 5.2: Sparsification Lemma (Lemma 6)
- Theorem 1: Approximate stabilizer rank bound
"""
function simulate_sum_over_cliffords(
    circuit::Vector{<:AbstractOperation}, 
    n_qubits::Int, 
    delta::Float64;
    verbose::Bool=false,
    rng::AbstractRNG=Random.GLOBAL_RNG)
    
    delta > 0 || throw(ArgumentError("Approximation error δ must be positive, got $delta"))
    n_qubits > 0 || throw(ArgumentError("Number of qubits must be positive, got $n_qubits"))
    
    non_clifford_count = count(op -> !isclifford(op), circuit)
    
    if non_clifford_count == 0
        if verbose
            @info "Pure Clifford circuit detected"
        end
        return simulate_pure_clifford_circuit(circuit, n_qubits, delta)
    end
    
    delta_per_gate = delta / non_clifford_count
    
    if verbose
        @info "Incremental Sum-over-Cliffords simulation"
        @info "  Non-Clifford gates: $non_clifford_count"
        @info "  Per-gate error budget: δᵢ = $(round(delta_per_gate, digits=6))"
    end
    
    initial_state = create_computational_zero_state(n_qubits)
    current_states = [initial_state]
    current_coefficients = [ComplexF64(1.0)]
    
    total_extent = 1.0
    non_clifford_idx = 0
    max_terms_seen = 1
    
    for (op_idx, op) in enumerate(circuit)
        if isclifford(op)
            for state in current_states
                apply!(state, op)
            end
        else
            non_clifford_idx += 1
            gate_decomp = get_gate_decomposition(op)
            gate_extent = stabilizer_extent(op)
            total_extent *= gate_extent
            num_decomp_terms = length(gate_decomp.coefficients)
            
            current_k = length(current_states)
            
            if verbose
                @info "Processing non-Clifford gate $non_clifford_idx/$non_clifford_count"
                @info "  Gate: $(typeof(op)), extent ξ = $(round(gate_extent, digits=4))"
                @info "  Decomposition terms: $num_decomp_terms"
                @info "  Current states before expansion: $current_k"
            end
            
            expanded_states = Vector{MixedDestabilizer}()
            expanded_coefficients = Vector{ComplexF64}()
            
            expected_size = current_k * num_decomp_terms
            sizehint!(expanded_states, expected_size)
            sizehint!(expanded_coefficients, expected_size)
            
            for (state, coeff) in zip(current_states, current_coefficients)
                for (gate_coeff, clifford_ops) in zip(gate_decomp.coefficients,
                                                       gate_decomp.clifford_operations)
                    new_state = copy(state)
                    for cliff_op in clifford_ops
                        apply!(new_state, cliff_op)
                    end
                    push!(expanded_states, new_state)
                    push!(expanded_coefficients, coeff * gate_coeff)
                end
            end
            
            expanded_count = length(expanded_states)
            max_terms_seen = max(max_terms_seen, expanded_count)
            
            if verbose
                @info "  States after expansion: $expanded_count"
                l1_before = sum(abs, expanded_coefficients)
                @info "  L1 norm: $(round(l1_before, digits=4))"
            end
            
            sparse_result = sparsify_mixed_destabilizer_decomposition(
                expanded_coefficients, expanded_states, delta_per_gate; rng=rng)
            
            current_states = sparse_result.states
            current_coefficients = sparse_result.coefficients
            
            if verbose
                @info "  States after sparsification: $(sparse_result.k)"
                @info "  Sparsification reduced by factor: $(round(expanded_count / sparse_result.k, digits=2))x"
            end
        end
    end
    
    if verbose
        @info "Simulation complete"
        @info "  Final number of states: $(length(current_states))"
        @info "  Maximum states during simulation: $max_terms_seen"
        @info "  Total extent ∏ξᵢ: $(round(total_extent, digits=4))"
        @info "  Expected final k bound: $(ceil(Int, 1 + total_extent/delta^2))"
    end
    
    final_stabilizer_states = Vector{Stabilizer}(undef, length(current_states))
    for (i, state) in enumerate(current_states)
        final_stabilizer_states[i] = Stabilizer(stabilizerview(state))
    end
    
    return SimulationState(
        final_stabilizer_states,
        current_coefficients,
        length(final_stabilizer_states),
        delta,
        total_extent
    )
end

"""
    simulate_pure_clifford_circuit(circuit, n_qubits, delta)

Handle pure Clifford circuits exactly (no approximation needed).
Returns a trivial SimulationState with single stabilizer state.
"""
function simulate_pure_clifford_circuit(circuit::AbstractVector{<:AbstractOperation}, 
                                       n_qubits::Int, 
                                       delta::Float64)
    state = create_computational_zero_state(n_qubits)
    
    for op in circuit
        try
            apply!(state, op)
        catch e
            @error "Failed to apply Clifford operation $op: $e"
            rethrow(e)
        end
    end
    
    if rank(state) != n_qubits
        throw(ErrorException(
            "Circuit produced mixed state (rank=$(rank(state))) from pure input - not a valid Clifford circuit"))
    end
    
    final_state = Stabilizer(stabilizerview(state))
    
    if size(final_state, 1) != n_qubits
        throw(ErrorException(
            "Internal error: stabilizer has $(size(final_state, 1)) generators, expected $n_qubits"))
    end
    
    return SimulationState([final_state], [ComplexF64(1.0)], 1, 0.0, 1.0)
end

"""
    MetropolisSamplerState

Cached state for efficient Metropolis sampling.
Stores current bitstring and precomputed amplitude to avoid recomputation.
"""
mutable struct MetropolisSamplerState
    current_bitstring::BitVector
    current_amplitude_sq::Float64
    states::Vector{Stabilizer}
    coefficients::Vector{ComplexF64}
    n_qubits::Int
end

"""
    initialize_metropolis(states, coeffs) -> MetropolisSamplerState

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
    metropolis_step!(metro::MetropolisSamplerState) -> Bool

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
    metropolis_mixing!(metro::MetropolisSamplerState, burn_in::Int)

Run burn-in period to reach equilibrium distribution.
Recommended: burn_in ≈ 10n for shallow circuits, 100n for deep circuits.
"""
function metropolis_mixing!(metro::MetropolisSamplerState, burn_in::Int)
    for _ in 1:burn_in
        metropolis_step!(metro)
    end
end

"""
    sample_measurement_outcomes(sim_state, n_samples; burn_in=0, verbose=false) -> Matrix{Bool}

MCMC Metropolis sampling from Section 4.2.
Complexity: O(k × n³ × T) where T = burn_in + n_samples × thinning.

# Arguments
- `sim_state::SimulationState`: Sparse stabilizer decomposition from sum-over-Cliffords
- `n_samples::Int`: Number of measurement outcomes to generate
- `burn_in::Int=0`: Mixing time (0 = auto-tune based on circuit)
- `verbose::Bool=false`: Show progress information

# Returns
`Matrix{Bool}` of size (n_samples, n_qubits) with measurement outcomes.
"""
function sample_measurement_outcomes(sim_state::SimulationState, 
                                    n_samples::Int; 
                                    burn_in::Int=0,
                                    verbose::Bool=false)
    
    n_samples > 0 || throw(ArgumentError("Number of samples must be positive, got $n_samples"))
    length(sim_state.sparse_states) > 0 || throw(ArgumentError("Simulation state must contain at least one stabilizer state"))
    
    n_qubits = nqubits(sim_state.sparse_states[1])
    
    if burn_in == 0
        circuit_depth_proxy = log2(max(1, sim_state.simulation_cost))
        burn_in = ceil(Int, 10 * n_qubits * (1 + circuit_depth_proxy / 10))
        if verbose
            @info "Auto-tuned burn-in: $burn_in steps"
        end
    end
    
    if verbose
        @info "Initializing Metropolis sampler with $(length(sim_state.sparse_states)) stabilizer states..."
    end
    metro = initialize_metropolis(sim_state.sparse_states, sim_state.coefficients)
    
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
    compute_outcome_frequencies(measurements::Matrix{Bool}) -> Dict{BitVector, Float64}

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
    validate_simulation_parameters(circuit, n_qubits, trajectories, delta)

 validation of simulation parameters.
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
    validate_operation_qubits(op, n_qubits, idx)

Validate that operation uses valid qubit indices.
"""
function validate_operation_qubits(op::AbstractOperation, n_qubits::Int, idx::Int)
    if op isa TGate
        1 ≤ op.qubit ≤ n_qubits || throw(ArgumentError(
            "Operation $idx: TGate qubit index $(op.qubit) out of range [1,$n_qubits]"))
        return nothing
    elseif op isa CCZGate
        for q in op.qubits
            1 ≤ q ≤ n_qubits || throw(ArgumentError(
                "Operation $idx: CCZGate qubit index $q out of range [1,$n_qubits]"))
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
    infer_circuit_nqubits(circuit) -> Int

Infer number of qubits from circuit operations.
"""
function infer_circuit_nqubits(circuit)
    max_qubit = 1
    for op in circuit
        if op isa TGate
            max_qubit = max(max_qubit, op.qubit)
        elseif op isa CCZGate
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
    lrtrajectories(circuit; trajectories=1000, delta=0.1, verbose=false) -> LRTrajectoryResults

Simulate quantum circuit with non-Clifford gates using low-rank stabilizer decomposition.

This is the main entry point for non-Clifford simulation, analogous to `pftrajectories` 
for Pauli frame simulation.

# Arguments
- `circuit`: Vector of quantum operations (Clifford gates, TGate, CCZGate, etc.)
- `trajectories::Int=1000`: Number of measurement samples to generate
- `delta::Float64=0.1`: Approximation error parameter (smaller = more accurate but slower)
- `verbose::Bool=false`: Show detailed progress information

# Returns
`LRTrajectoryResults` containing measurement outcomes and simulation statistics.
Use `lrmeasurements(result)` to extract the measurement matrix.

# Algorithm
Implements the Sum-over-Cliffords method from Section 2.3.2 of Bravyi et al. 2019
with incremental sparsification (Section 5.2) after each non-Clifford gate:

1. Initialize with |0⟩ state
2. For each gate:
   - Clifford: apply to all current stabilizer states
   - Non-Clifford: expand decomposition, then sparsify immediately
3. Sample measurement outcomes using Metropolis MCMC (Section 4.2)

# Performance
- Memory: O(k × n) where k ≈ total_extent/δ² (never stores 2^m terms)
- Simulation cost scales as ∏ⱼ ξ(Vⱼ) where ξ(T) ≈ 1.172, ξ(CCZ) = 16/9
- Time complexity: O(m × k² × n³) for simulation + O(k × n³ × trajectories) for sampling

# Example
```julia
# Simple T-gate circuit
circuit = [sHadamard(1), TGate(1), sHadamard(1)]
result = lrtrajectories(circuit; trajectories=5000, delta=0.1)
println(result)  # Shows statistics and top outcomes
measurements = lrmeasurements(result)  # Extract raw data

# Check probability distribution
p0 = sum(measurements[:, 1] .== false) / size(measurements, 1)
println("P(|0⟩) = \$p0")  # Should be ≈ cos²(π/8) ≈ 0.854
```
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
    
    if verbose
        @info "Starting low-rank stabilizer simulation"
        @info "  Qubits: $n_qubits"
        @info "  Gates: $(length(circuit))"
        @info "  Trajectories: $trajectories"
        @info "  δ: $delta"
    end
    
    circuit_vec = collect(AbstractOperation, circuit)
    
    try
        if verbose
            @info "Step 1/2: Computing sum-over-Cliffords decomposition..."
        end
        sim_state = simulate_sum_over_cliffords(circuit_vec, n_qubits, delta; verbose)
        
        if verbose
            @info "Decomposition complete: $(sim_state.simulation_cost) sparse terms"
            @info "Total extent: $(sim_state.original_extent)"
        end
        
        if verbose
            @info "Step 2/2: Sampling measurement outcomes..."
        end
        measurements = sample_measurement_outcomes(sim_state, trajectories; verbose)
        
        total_runtime = time() - start_time
        
        if verbose
            @info "Simulation completed in $(round(total_runtime, digits=2)) seconds"
        end
        
        return LRTrajectoryResults(
            measurements,
            sim_state.simulation_cost,
            sim_state.approximation_error,
            sim_state.original_extent,
            n_qubits,
            total_runtime
        )
        
    catch e
        @error "Simulation failed: $e"
        rethrow(e)
    end
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
    lrcost(circuit; delta=0.1) -> NamedTuple

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
```julia
circuit = [sHadamard(1), TGate(1), TGate(1), TGate(1), sHadamard(1)]
cost = lrcost(circuit; delta=0.1)
println("Estimated cost: \$(cost.estimated_k) stabilizer terms")
println("Total extent: \$(cost.total_extent)")

if cost.estimated_k > 100000
    println("Consider using larger delta for faster simulation")
end
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