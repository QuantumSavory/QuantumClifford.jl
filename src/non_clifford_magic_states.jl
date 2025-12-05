# Implements Lemma 6 (Sparsification), Lemma 1 (Lifting), and Proposition 2

using Random
using Statistics

"""
    SparsifiedState

Result of applying Sparsification Lemma (Lemma 6) from Section 5.2.
Represents |Ω⟩ = (||c||₁/k) Σₐ₌₁ᵏ |ωₐ⟩ that approximates |ψ⟩ = Σⱼ cⱼ|φⱼ⟩.
"""
struct SparsifiedState
    states::Vector{<:Stabilizer}
    coefficients::Vector{ComplexF64}
    k::Int
    original_l1_norm::Float64
    approximation_error_bound::Float64
    expected_norm_bound::Float64
    
    function SparsifiedState(states, coeffs, k, l1_norm, delta)
        length(states) == length(coeffs) == k || throw(DimensionMismatch("Inconsistent array lengths"))
        expected_norm = 1.0 + l1_norm^2 / k
        new(states, coeffs, k, l1_norm, delta, expected_norm)
    end
end

"""
    sparsify_stabilizer_decomposition(coefficients, states, delta; rng=Random.GLOBAL_RNG)

Implements Sparsification Lemma (Lemma 6) from Section 5.2.

Given |ψ⟩ = Σⱼ cⱼ|φⱼ⟩ with ||c||₁, constructs random state |Ω⟩ = (||c||₁/k) Σₐ₌₁ᵏ |ωₐ⟩
where each |ωₐ⟩ is sampled from {|φⱼ⟩} with probability |cⱼ|/||c||₁.

Returns approximation satisfying E[||ψ - Ω||²] = ||c||₁²/k ≤ δ² for k ≥ ||c||₁²/δ².
"""
function sparsify_stabilizer_decomposition(coefficients::Vector{ComplexF64}, 
                                         states::Vector{<:Stabilizer}, 
                                         delta::Float64;
                                         rng::AbstractRNG=Random.GLOBAL_RNG)
    
    length(coefficients) == length(states) || throw(DimensionMismatch("Coefficients and states must have same length"))
    delta > 0 || throw(ArgumentError("Approximation error δ must be positive"))
    
    l1_norm = sum(abs.(coefficients))
    l1_norm > 0 || throw(ArgumentError("All coefficients are zero - no valid decomposition"))
    
    k = max(1, ceil(Int, l1_norm^2 / delta^2))
    
    abs_coeffs = abs.(coefficients)
    probabilities = abs_coeffs / l1_norm
    cumulative_probs = cumsum(probabilities)
    
    sampled_states = typeof(states)()
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
    
    return SparsifiedState(sampled_states, sampled_coefficients, k, l1_norm, delta)
end

"""
    estimate_sparsification_quality(sparse::SparsifiedState)

Estimate quality bounds from Lemma 7 (Sparsification tail bound).
"""
function estimate_sparsification_quality(sparse::SparsifiedState)
    return (
        k=sparse.k,
        expected_error=sparse.original_l1_norm^2 / sparse.k,
        error_bound=sparse.approximation_error_bound,
        expected_norm=sparse.expected_norm_bound
    )
end

"""
    MagicStateDecomposition

Stabilizer decomposition V|+⟩^⊗t = Σⱼ cⱼ|φⱼ⟩ for Clifford magic states.
Used as intermediate step for Lifting Lemma.

For Clifford magic states: ξ(ψ) = F(ψ)⁻¹ where F(ψ) = max_φ |⟨φ|ψ⟩|² (Proposition 2).
"""
struct MagicStateDecomposition
    gate_type::Symbol
    coefficients::Vector{ComplexF64}
    stabilizer_states::Vector{<:Stabilizer}
    l1_norm::Float64
    stabilizer_extent::Float64
    stabilizer_fidelity::Float64
    
    function MagicStateDecomposition(gate_type, coeffs, states)
        length(coeffs) == length(states) || throw(DimensionMismatch("Mismatched coefficients and states"))
        l1 = sum(abs.(coeffs))
        
        if gate_type == :CCZ
            xi = 16.0/9.0
            fidelity = 9.0/16.0
        elseif gate_type == :R_theta
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
    CliffordGateDecomposition

Sum-over-Cliffords decomposition U = Σⱼ cⱼKⱼ where Kⱼ are Clifford unitaries.
"""
struct CliffordGateDecomposition
    gate_type::Symbol
    coefficients::Vector{ComplexF64}
    clifford_operations::Vector{Vector{<:AbstractOperation}}
    l1_norm::Float64
    stabilizer_extent::Float64
    target_qubits::Vector{Int}
    
    function CliffordGateDecomposition(gate_type, coeffs, ops, qubits)
        length(coeffs) == length(ops) || throw(DimensionMismatch("Mismatched coefficients and operations"))
        l1 = sum(abs.(coeffs))
        xi = l1^2
        new(gate_type, coeffs, ops, l1, xi, qubits)
    end
end

"""
    decompose_rotation_magic_state(θ; nqubits=1)

Create magic state decomposition R(θ)|+⟩ = Σⱼ cⱼ|φⱼ⟩ from Eq. (26).

R(θ)|+⟩ = (cos(θ/2) - sin(θ/2))|+⟩ + √2 sin(θ/2)e^(-iπ/4)S|+⟩

Returns optimal decomposition with ξ(R(θ)) = (cos(θ/2) + tan(π/8)sin(θ/2))².
"""
function decompose_rotation_magic_state(θ::Float64; nqubits::Int=1)
    cos_half = cos(θ/2)
    sin_half = sin(θ/2)
    
    c1 = ComplexF64(cos_half - sin_half, 0.0)
    c2 = sqrt(2) * sin_half * exp(-im * π/4)
    
    if nqubits == 1
        plus_state = Stabilizer([P"X"])
        s_plus_state = Stabilizer([P"Y"])
    else
        plus_generators = [PauliOperator(nqubits) for _ in 1:nqubits]
        for i in 1:nqubits
            plus_generators[i] = zero(PauliOperator, nqubits)
            plus_generators[i][i] = (true, false)
        end
        plus_state = Stabilizer(plus_generators)
        
        s_plus_generators = [PauliOperator(nqubits) for _ in 1:nqubits]
        for i in 1:nqubits
            s_plus_generators[i] = zero(PauliOperator, nqubits)
            s_plus_generators[i][i] = (true, true)
        end
        s_plus_state = Stabilizer(s_plus_generators)
    end
    
    return MagicStateDecomposition(:R_theta, [c1, c2], [plus_state, s_plus_state])
end

"""
    lifting_lemma_single_qubit(magic_decomp::MagicStateDecomposition, qubit::Int)

Apply Lifting Lemma (Lemma 1) to convert R(θ)|+⟩ = Σⱼ cⱼ|φⱼ⟩ to R(θ) = Σⱼ cⱼKⱼ.

For diagonal gate, if V|+⟩ = Σⱼ cⱼKⱼ|+⟩ where Kⱼ are diagonal Cliffords, then V = Σⱼ cⱼKⱼ.
"""
function lifting_lemma_single_qubit(magic_decomp::MagicStateDecomposition, qubit::Int)
    coeffs = magic_decomp.coefficients
    
    clifford_ops = Vector{Vector{AbstractOperation}}()
    
    if length(coeffs) == 2
        push!(clifford_ops, AbstractOperation[])
        push!(clifford_ops, [sPhase(qubit)])
    else
        throw(ArgumentError("Unexpected number of terms in single-qubit decomposition"))
    end
    
    return CliffordGateDecomposition(magic_decomp.gate_type, coeffs, clifford_ops, [qubit])
end

"""
    decompose_T_gate(qubit::Int)

T gate decomposition: T = R(π/4) with ξ(T) = (cos(π/8) + tan(π/8)sin(π/8))².
"""
function decompose_T_gate(qubit::Int)
    magic_decomp = decompose_rotation_magic_state(π/4; nqubits=1)
    gate_decomp = lifting_lemma_single_qubit(magic_decomp, qubit)
    return CliffordGateDecomposition(:T, gate_decomp.coefficients, 
                                   gate_decomp.clifford_operations, [qubit])
end

"""
    create_ccz_stabilizer_state(operations::Vector{<:AbstractOperation})

Create stabilizer state by applying Clifford operations to |+++⟩.
This computes the effect of CZ and Z operations on the stabilizer generators.
"""
function create_ccz_stabilizer_state(operations::Vector{<:AbstractOperation})
    state = MixedDestabilizer(S"XII IXI IIX")
    
    for op in operations
        apply!(state, op)
    end
    
    return Stabilizer(stabilizerview(state))
end

"""
    decompose_CCZ_magic_state()

Create optimal magic state decomposition for CCZ gate using Proposition 2.

For Clifford magic states: ξ(ψ) = F(ψ)⁻¹ where F(ψ) = max_φ |⟨φ|ψ⟩|².
For CCZ: F(CCZ) = |⟨+++|CCZ⟩|² = 9/16, so ξ(CCZ) = 16/9.

Uses group decomposition |CCZ⟩ = (1/|Q|⟨CCZ|+++⟩) Σ_{q∈Q} q|+++⟩
where Q = ⟨X₁CZ₂,₃, X₂CZ₁,₃, X₃CZ₁,₂⟩ has 8 elements.

Returns optimal decomposition with ξ(CCZ) = 16/9.
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
    
    return MagicStateDecomposition(:CCZ, coefficients, states)
end

"""
    lifting_lemma_CCZ(magic_decomp::MagicStateDecomposition, qubits::Vector{Int})

Apply Lifting Lemma to convert CCZ magic state decomposition to gate decomposition.
"""
function lifting_lemma_CCZ(magic_decomp::MagicStateDecomposition, qubits::Vector{Int})
    length(qubits) == 3 || throw(ArgumentError("CCZ requires exactly 3 qubits"))
    
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
    
    return CliffordGateDecomposition(:CCZ, coeffs, clifford_ops, qubits)
end

"""
    decompose_CCZ_gate(qubits::Vector{Int})

Get optimal CCZ gate decomposition with ξ(CCZ) = 16/9
"""
function decompose_CCZ_gate(qubits::Vector{Int})
    magic_decomp = decompose_CCZ_magic_state()
    return lifting_lemma_CCZ(magic_decomp, qubits)
end