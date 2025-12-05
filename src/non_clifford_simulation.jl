# Sum-over-Cliffords circuit simulation for non-Clifford gates
# Implements simulation methods from Section 2.3

"""
    TGate <: AbstractOperation

Marker for T-gate (π/8 phase rotation).
Triggers optimal decomposition with ξ(T) ≈ 1.17 in simulation.
"""
struct TGate <: AbstractOperation
    qubit::Int
    
    function TGate(q::Int)
        q > 0 || throw(ArgumentError("Qubit index must be positive"))
        new(q)
    end
end

"""
    CCZGate <: AbstractOperation

Marker for Controlled-Controlled-Z gate.
Triggers optimal decomposition with ξ(CCZ) = 16/9 ≈ 1.78 in simulation.
"""
struct CCZGate <: AbstractOperation
    qubits::Vector{Int}
    
    function CCZGate(q1::Int, q2::Int, q3::Int)
        all(q -> q > 0, [q1, q2, q3]) || throw(ArgumentError("All qubit indices must be positive"))
        length(unique([q1, q2, q3])) == 3 || throw(ArgumentError("CCZ requires 3 distinct qubits"))
        new(sort([q1, q2, q3]))
    end
end

nqubits(::TGate) = 1
nqubits(::CCZGate) = 3

"""
    identify_gate_type(op::AbstractOperation)

Gate classification using direct isa checks.
Identifies all Clifford operations and distinguishes parametric non-Clifford gates.
"""
function identify_gate_type(op::AbstractOperation)
    if op isa TGate || op isa CCZGate
        return :non_clifford
    end
    
    if op isa sHadamard || op isa sPhase || op isa sInvPhase || 
       op isa sX || op isa sY || op isa sZ || op isa sId1 ||
       op isa sSQRTX || op isa sInvSQRTX || op isa sSQRTY || op isa sInvSQRTY ||
       op isa sHadamardXY || op isa sHadamardYZ || op isa sCXYZ || op isa sCZYX
        return :clifford
    end
    
    if op isa sCNOT || op isa sCPHASE || op isa sSWAP ||
       op isa sXCX || op isa sXCY || op isa sXCZ ||
       op isa sYCX || op isa sYCY || op isa sYCZ ||
       op isa sZCX || op isa sZCY || op isa sZCZ ||
       op isa sSWAPCX || op isa sInvSWAPCX || op isa sCZSWAP || op isa sCXSWAP ||
       op isa sISWAP || op isa sInvISWAP || op isa sSQRTZZ || op isa sInvSQRTZZ
        return :clifford
    end
    
    if op isa AbstractSingleQubitOperator || op isa AbstractTwoQubitOperator
        return :clifford
    end
    
    if op isa AbstractCliffordOperator || op isa AbstractSymbolicOperator
        return :clifford
    end
    
    if op isa sMX || op isa sMY || op isa sMZ || op isa sMRX || op isa sMRY || op isa sMRZ
        return :clifford
    end
    
    if op isa PauliMeasurement
        return :clifford
    end
    
    if op isa SparseGate
        return identify_gate_type(op.cliff)
    end
    
    return :non_clifford
end

"""
    extract_gate_parameters(op::AbstractOperation)

Extract parameters from gates using QuantumClifford.jl's structure.
Returns (gate_type, parameters, qubits) for both Clifford and non-Clifford gates.
"""
function extract_gate_parameters(op::AbstractOperation)
    if op isa TGate
        return (:T, [π/4], [op.qubit])
    elseif op isa CCZGate
        return (:CCZ, Float64[], op.qubits)
    end
    
    op_type = typeof(op)
    
    if hasfield(op_type, :q) && !hasfield(op_type, :q2)
        qubit = op.q
        qubits = [qubit]
        
        if op isa sHadamard
            return (:hadamard, Float64[], qubits)
        elseif op isa sPhase
            return (:phase, [π/2], qubits)
        elseif op isa sInvPhase  
            return (:inv_phase, [-π/2], qubits)
        elseif op isa sX
            return (:pauli_x, Float64[], qubits)
        elseif op isa sY
            return (:pauli_y, Float64[], qubits)
        elseif op isa sZ
            return (:pauli_z, Float64[], qubits)
        elseif op isa sId1
            return (:identity, Float64[], qubits)
        elseif op isa sSQRTX
            return (:sqrt_x, [π/2], qubits)
        elseif op isa sInvSQRTX
            return (:inv_sqrt_x, [-π/2], qubits)
        elseif op isa sSQRTY
            return (:sqrt_y, [π/2], qubits)
        elseif op isa sInvSQRTY
            return (:inv_sqrt_y, [-π/2], qubits)
        else
            return (:T, [π/4], qubits)
        end
    
    elseif hasfield(op_type, :q1) && hasfield(op_type, :q2)
        qubits = [op.q1, op.q2]
        
        if op isa sCNOT
            return (:cnot, Float64[], qubits)
        elseif op isa sCPHASE
            return (:cphase, [π], qubits)
        elseif op isa sSWAP
            return (:swap, Float64[], qubits)
        elseif op isa sXCX || op isa sXCY || op isa sXCZ ||
               op isa sYCX || op isa sYCY || op isa sYCZ ||
               op isa sZCX || op isa sZCY || op isa sZCZ
            return (:controlled_pauli, Float64[], qubits)
        elseif op isa sISWAP
            return (:iswap, [π/2], qubits)
        elseif op isa sInvISWAP
            return (:inv_iswap, [-π/2], qubits)
        else
            return (:unknown_two_qubit, Float64[], qubits)
        end
    
    elseif op isa sMX || op isa sMY || op isa sMZ
        qubit = hasfield(typeof(op), :qubit) ? op.qubit : 1
        measurement_basis = op isa sMX ? :x : (op isa sMY ? :y : :z)
        return (Symbol(:measure_, measurement_basis), Float64[], [qubit])
    
    elseif op isa SparseGate
        gate_type, params, _ = extract_gate_parameters(op.cliff)
        return (gate_type, params, op.indices)
    
    elseif op isa PauliMeasurement
        n_qubits = nqubits(op.pauli)
        return (:pauli_measurement, Float64[], collect(1:n_qubits))
    
    else
        qubits = if hasfield(op_type, :qubits)
            op.qubits
        elseif hasfield(op_type, :indices) 
            op.indices
        else
            [1]
        end
        
        return (:unknown, Float64[], qubits)
    end
end

"""
    get_optimal_gate_decomposition(gate_type::Symbol, parameters, qubits)

Get optimal CliffordGateDecomposition for all supported gate types.
"""
function get_optimal_gate_decomposition(gate_type::Symbol, parameters, qubits)
    if gate_type == :T
        return decompose_T_gate(qubits[1])
    elseif gate_type == :phase && length(parameters) == 1
        θ = parameters[1]
        magic_decomp = decompose_rotation_magic_state(θ)
        return lifting_lemma_single_qubit(magic_decomp, qubits[1])
    elseif gate_type == :CCZ
        length(qubits) == 3 || throw(ArgumentError("CCZ requires exactly 3 qubits"))
        return decompose_CCZ_gate(qubits)
    elseif gate_type in [:sqrt_x, :inv_sqrt_x, :sqrt_y, :inv_sqrt_y] && length(parameters) == 1
        θ = parameters[1]
        return CliffordGateDecomposition(gate_type, [ComplexF64(1.0)], 
                                       [AbstractOperation[]], qubits)
    elseif gate_type in [:hadamard, :pauli_x, :pauli_y, :pauli_z, :identity]
        return CliffordGateDecomposition(gate_type, [ComplexF64(1.0)], 
                                       [AbstractOperation[]], qubits)
    elseif gate_type in [:cnot, :cphase, :swap, :controlled_pauli]
        return CliffordGateDecomposition(gate_type, [ComplexF64(1.0)], 
                                       [AbstractOperation[]], qubits)
    elseif gate_type == :unknown_two_qubit && length(qubits) == 3
        return decompose_CCZ_gate(qubits)
    else
        @warn "Unknown gate type $gate_type, using T-gate decomposition as fallback"
        return decompose_T_gate(qubits[1])
    end
end

"""
    CircuitDecomposition

Full sum-over-Cliffords representation U = Σⱼ cⱼKⱼ where Kⱼ are Clifford circuits.
"""
struct CircuitDecomposition
    coefficients::Vector{ComplexF64}
    clifford_circuits::Vector{Vector{<:AbstractOperation}}
    l1_norm::Float64
    stabilizer_extent::Float64
    n_qubits::Int
    
    function CircuitDecomposition(coeffs, circuits, n_qubits)
        length(coeffs) == length(circuits) || throw(DimensionMismatch("Mismatched coefficients and circuits"))
        l1 = sum(abs.(coeffs))
        xi = l1^2
        new(coeffs, circuits, l1, xi, n_qubits)
    end
end

"""
    SimulationResult

Final sparse representation ready for sampling and probability estimation.
"""
struct SimulationResult
    sparse_states::Vector{<:Stabilizer}
    coefficients::Vector{ComplexF64}
    simulation_cost::Int
    approximation_error::Float64
    original_extent::Float64
end

"""
    create_computational_zero_state(n_qubits::Int)

Creates proper |0ⁿ⟩ state stabilized by Z₁, Z₂, ..., Zₙ.
Used as initial state for Sum-over-Cliffords method.
"""
function create_computational_zero_state(n_qubits::Int)::MixedDestabilizer
    n_qubits > 0 || throw(ArgumentError("Number of qubits must be positive"))
    
    z_operators = [zero(PauliOperator, n_qubits) for _ in 1:n_qubits]
    for i in 1:n_qubits
        z_operators[i][i] = (false, true)
    end
    
    return MixedDestabilizer(Stabilizer(z_operators))
end

"""
    construct_full_circuit_decomposition(clifford_sections, gate_decompositions)

Full implementation of circuit combination from Section 2.3.2.
"""
function construct_full_circuit_decomposition(clifford_sections::Vector{Vector{AbstractOperation}}, 
                                            gate_decompositions::Vector{CliffordGateDecomposition})
    final_coeffs = ComplexF64[]
    final_circuits = Vector{Vector{AbstractOperation}}()
    
    num_gates = length(gate_decompositions)
    if num_gates == 0
        push!(final_coeffs, ComplexF64(1.0))
        push!(final_circuits, vcat(clifford_sections...))
        return final_coeffs, final_circuits
    end
    
    decomp_sizes = [length(decomp.coefficients) for decomp in gate_decompositions]
    
    for indices in Iterators.product([1:size for size in decomp_sizes]...)
        coeff = ComplexF64(1.0)
        for (j, k) in enumerate(indices)
            coeff *= gate_decompositions[j].coefficients[k]
        end
        
        circuit = AbstractOperation[]
        append!(circuit, clifford_sections[1])
        
        for (j, k) in enumerate(indices)
            append!(circuit, gate_decompositions[j].clifford_operations[k])
            append!(circuit, clifford_sections[j+1])
        end
        
        push!(final_coeffs, coeff)
        push!(final_circuits, circuit)
    end
    
    return final_coeffs, final_circuits
end

"""
    simulate_sum_over_cliffords(circuit, n_qubits, delta)

Full implementation of Sum-over-Cliffords simulation method from Section 2.3.2.
"""
function simulate_sum_over_cliffords(circuit::Vector{<:AbstractOperation}, 
                                   n_qubits::Int, 
                                   delta::Float64)
    delta > 0 || throw(ArgumentError("Approximation error δ must be positive"))
    n_qubits > 0 || throw(ArgumentError("Number of qubits must be positive"))
    
    clifford_sections = Vector{Vector{AbstractOperation}}()
    non_clifford_gates = Vector{AbstractOperation}()
    
    current_clifford_section = AbstractOperation[]
    
    for op in circuit
        if identify_gate_type(op) == :clifford
            push!(current_clifford_section, op)
        else
            push!(clifford_sections, copy(current_clifford_section))
            push!(non_clifford_gates, op)
            current_clifford_section = AbstractOperation[]
        end
    end
    push!(clifford_sections, current_clifford_section)
    
    if isempty(non_clifford_gates)
        return simulate_pure_clifford_circuit(circuit, n_qubits, delta)
    end
    
    gate_decompositions = CliffordGateDecomposition[]
    
    for gate in non_clifford_gates
        gate_type, parameters, qubits = extract_gate_parameters(gate)
        decomp = get_optimal_gate_decomposition(gate_type, parameters, qubits)
        push!(gate_decompositions, decomp)
    end
    
    final_coeffs, final_circuits = construct_full_circuit_decomposition(
        clifford_sections, gate_decompositions)
    
    stabilizer_states = Vector{Stabilizer}()
    
    for clifford_circuit in final_circuits
        state = create_computational_zero_state(n_qubits)
        
        for op in clifford_circuit
            apply!(state, op)
        end
        
        push!(stabilizer_states, Stabilizer(stabilizerview(state)))
    end
    
    sparse_result = sparsify_stabilizer_decomposition(
        final_coeffs, 
        stabilizer_states, 
        delta
    )
    
    total_l1_norm = sum(abs.(final_coeffs))
    total_extent = total_l1_norm^2
    
    return SimulationResult(
        sparse_result.states,
        sparse_result.coefficients,
        sparse_result.k,
        delta,
        total_extent
    )
end

"""
    simulate_pure_clifford_circuit(circuit, n_qubits, delta)

Handle pure Clifford circuits exactly (no approximation needed).
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
        throw(ErrorException("Circuit produced mixed state (rank=$(rank(state))) from pure input - not a valid Clifford circuit"))
    end
    
    final_state = Stabilizer(stabilizerview(state))
    
    if size(final_state, 1) != n_qubits
        throw(ErrorException("Internal error: stabilizer has $(size(final_state, 1)) generators, expected $n_qubits"))
    end
    
    return SimulationResult([final_state], [ComplexF64(1.0)], 1, 0.0, 1.0)
end

"""
    estimate_simulation_cost(circuit, delta)

Estimate simulation cost before running full simulation.
"""
function estimate_simulation_cost(circuit::Vector{<:AbstractOperation}, delta::Float64)
    total_extent = 1.0
    non_clifford_count = 0
    gate_extents = Float64[]
    
    for op in circuit
        if identify_gate_type(op) == :non_clifford
            gate_type, parameters, qubits = extract_gate_parameters(op)
            
            if gate_type == :T
                xi = (cos(π/8) + tan(π/8) * sin(π/8))^2
            elseif gate_type == :CCZ
                xi = 16.0/9.0
            elseif gate_type == :phase && length(parameters) == 1
                θ = parameters[1]
                xi = (cos(θ/2) + tan(π/8) * sin(θ/2))^2
            else
                xi = 2.0
            end
            
            total_extent *= xi
            non_clifford_count += 1
            push!(gate_extents, xi)
        end
    end
    
    k_bound = if non_clifford_count == 0
        1
    else
        ceil(Int, 1 + total_extent / delta^2)
    end
    
    return (
        estimated_cost=k_bound,
        total_extent=total_extent,
        non_clifford_gates=non_clifford_count,
        individual_extents=gate_extents,
        scaling_factor=non_clifford_count > 0 ? total_extent^(1/non_clifford_count) : 1.0
    )
end