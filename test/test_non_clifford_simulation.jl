@testitem "Non-Clifford Circuit Simulation" tags=[:non_clifford] begin
    using QuantumClifford
    import QuantumClifford: AbstractOperation, nqubits, identify_gate_type, extract_gate_parameters, CircuitDecomposition, get_optimal_gate_decomposition, estimate_simulation_cost

    @testset "Pure Clifford Circuit" begin
        circuit_clifford = [sHadamard(1), sCNOT(1,2), sZ(1)]
        
        result = simulate_non_clifford_circuit(circuit_clifford, 2, samples=50, delta=0.1, verbose=false)
        @test result.simulation_cost == 1
        @test result.approximation_error == 0.0
        @test length(result.outcomes) == 50
    end

    @testset "Mixed Circuit" begin
        circuit_mixed = [sHadamard(1), sPhase(1), sCNOT(1,2)]
        
        result = simulate_non_clifford_circuit(circuit_mixed, 2, samples=30, delta=0.2, verbose=false)
        @test result.simulation_cost >= 1
        @test length(result.outcomes) == 30
    end

    @testset "T-Gate Simulation" begin
        circuit_T = [sHadamard(1), TGate(1), sHadamard(1)]
        result_T = simulate_non_clifford_circuit(circuit_T, 1, samples=100, delta=0.1, verbose=false)
        
        expected_xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        expected_cost_T = ceil(Int, 1 + expected_xi_T / 0.1^2)
        
        @test result_T.simulation_cost > 1
        @test result_T.simulation_cost <= expected_cost_T
        @test result_T.approximation_error ≈ 0.1 atol=0.01
    end

    @testset "Multiple T-Gates" begin
        circuit_3T = [TGate(1), TGate(1), TGate(1)]
        result_3T = simulate_non_clifford_circuit(circuit_3T, 1, samples=50, delta=0.1, verbose=false)
        
        expected_xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        expected_xi_3T = expected_xi_T^3
        expected_cost_3T = ceil(Int, 1 + expected_xi_3T / 0.1^2)
        
        @test result_3T.simulation_cost <= expected_cost_3T
    end

    @testset "CCZ Gate Simulation" begin
        circuit_CCZ = [
            sHadamard(1), sHadamard(2), sHadamard(3),
            CCZGate(1, 2, 3),
            sHadamard(1), sHadamard(2), sHadamard(3)
        ]
        result_CCZ = simulate_non_clifford_circuit(circuit_CCZ, 3, samples=50, delta=0.1, verbose=false)
        
        expected_xi_CCZ = 16.0/9.0
        expected_cost_CCZ = ceil(Int, 1 + expected_xi_CCZ / 0.1^2)
        
        @test result_CCZ.simulation_cost > 1
        @test result_CCZ.simulation_cost <= expected_cost_CCZ
    end

    @testset "Error Handling" begin
        @test_throws ArgumentError simulate_non_clifford_circuit([], 1)
        @test_throws ArgumentError simulate_non_clifford_circuit([sHadamard(1)], 0)
        @test_throws ArgumentError simulate_non_clifford_circuit([sHadamard(1)], 1, samples=-1)
        @test_throws ArgumentError simulate_non_clifford_circuit([sHadamard(1)], 1, delta=1.5)
        
        circuit_invalid = [sHadamard(5)]
        @test_throws ArgumentError simulate_non_clifford_circuit(circuit_invalid, 2)
    end

    @testset "Cost Estimation" begin
        circuit = [sHadamard(1), TGate(1)]
        deltas = [0.5, 0.2, 0.1, 0.05]
        
        for delta in deltas
            result = simulate_non_clifford_circuit(circuit, 1, samples=30, delta=delta, verbose=false)
            xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
            bound = ceil(Int, 1 + xi_T / delta^2)
            
            @test result.simulation_cost <= bound
            @test result.approximation_error ≈ delta atol=0.02
        end
    end
    
    @testset "nqubits Interface" begin
        t = TGate(1)
        @test nqubits(t) == 1
        
        ccz = CCZGate(1, 2, 3)
        @test nqubits(ccz) == 3
    end
    
    @testset "Gate Type Identification - Single Qubit" begin
        @test identify_gate_type(TGate(1)) == :non_clifford
        
        @test identify_gate_type(sHadamard(1)) == :clifford
        @test identify_gate_type(sPhase(1)) == :clifford
        @test identify_gate_type(sInvPhase(1)) == :clifford
        @test identify_gate_type(sX(1)) == :clifford
        @test identify_gate_type(sY(1)) == :clifford
        @test identify_gate_type(sZ(1)) == :clifford
        @test identify_gate_type(sId1(1)) == :clifford
        @test identify_gate_type(sSQRTX(1)) == :clifford
        @test identify_gate_type(sInvSQRTX(1)) == :clifford
        @test identify_gate_type(sSQRTY(1)) == :clifford
        @test identify_gate_type(sInvSQRTY(1)) == :clifford
        @test identify_gate_type(sHadamardXY(1)) == :clifford
        @test identify_gate_type(sHadamardYZ(1)) == :clifford
        @test identify_gate_type(sCXYZ(1)) == :clifford
        @test identify_gate_type(sCZYX(1)) == :clifford
    end
    
    @testset "Gate Type Identification - Two Qubit" begin
        @test identify_gate_type(CCZGate(1, 2, 3)) == :non_clifford
        
        @test identify_gate_type(sCNOT(1, 2)) == :clifford
        @test identify_gate_type(sCPHASE(1, 2)) == :clifford
        @test identify_gate_type(sSWAP(1, 2)) == :clifford
        @test identify_gate_type(sXCX(1, 2)) == :clifford
        @test identify_gate_type(sXCY(1, 2)) == :clifford
        @test identify_gate_type(sXCZ(1, 2)) == :clifford
        @test identify_gate_type(sYCX(1, 2)) == :clifford
        @test identify_gate_type(sYCY(1, 2)) == :clifford
        @test identify_gate_type(sYCZ(1, 2)) == :clifford
        @test identify_gate_type(sZCX(1, 2)) == :clifford
        @test identify_gate_type(sZCY(1, 2)) == :clifford
        @test identify_gate_type(sZCZ(1, 2)) == :clifford
        @test identify_gate_type(sSWAPCX(1, 2)) == :clifford
        @test identify_gate_type(sInvSWAPCX(1, 2)) == :clifford
        @test identify_gate_type(sCZSWAP(1, 2)) == :clifford
        @test identify_gate_type(sCXSWAP(1, 2)) == :clifford
        @test identify_gate_type(sISWAP(1, 2)) == :clifford
        @test identify_gate_type(sInvISWAP(1, 2)) == :clifford
        @test identify_gate_type(sSQRTZZ(1, 2)) == :clifford
        @test identify_gate_type(sInvSQRTZZ(1, 2)) == :clifford
    end
    
    @testset "Gate Type Identification - Measurements" begin
        @test identify_gate_type(sMX(1, 1)) == :clifford
        @test identify_gate_type(sMY(1, 1)) == :clifford
        @test identify_gate_type(sMZ(1, 1)) == :clifford
        @test identify_gate_type(sMRX(1, 1)) == :clifford
        @test identify_gate_type(sMRY(1, 1)) == :clifford
        @test identify_gate_type(sMRZ(1, 1)) == :clifford
        @test identify_gate_type(PauliMeasurement(P"X", 1)) == :clifford
    end
    
    @testset "Gate Type Identification - SparseGate" begin
        sparse_clifford = SparseGate(tCNOT, [1, 2])
        @test identify_gate_type(sparse_clifford) == :clifford
    end

    
    @testset "Extract T and CCZ Parameters" begin
        gate_type, params, qubits = extract_gate_parameters(TGate(3))
        @test gate_type == :T
        @test params ≈ [π/4]
        @test qubits == [3]
        
        gate_type, params, qubits = extract_gate_parameters(CCZGate(1, 2, 3))
        @test gate_type == :CCZ
        @test params == Float64[]
        @test qubits == [1, 2, 3]
    end
    
    @testset "Extract Single Qubit Clifford Parameters" begin
        @test extract_gate_parameters(sHadamard(2)) == (:hadamard, Float64[], [2])
        @test extract_gate_parameters(sPhase(2))[1:2] == (:phase, [π/2])
        @test extract_gate_parameters(sInvPhase(2))[1:2] == (:inv_phase, [-π/2])
        @test extract_gate_parameters(sX(2)) == (:pauli_x, Float64[], [2])
        @test extract_gate_parameters(sY(2)) == (:pauli_y, Float64[], [2])
        @test extract_gate_parameters(sZ(2)) == (:pauli_z, Float64[], [2])
        @test extract_gate_parameters(sId1(2)) == (:identity, Float64[], [2])
        @test extract_gate_parameters(sSQRTX(2))[1:2] == (:sqrt_x, [π/2])
        @test extract_gate_parameters(sInvSQRTX(2))[1:2] == (:inv_sqrt_x, [-π/2])
        @test extract_gate_parameters(sSQRTY(2))[1:2] == (:sqrt_y, [π/2])
        @test extract_gate_parameters(sInvSQRTY(2))[1:2] == (:inv_sqrt_y, [-π/2])
    end
    
    @testset "Extract Two Qubit Clifford Parameters" begin
        @test extract_gate_parameters(sCNOT(1, 2)) == (:cnot, Float64[], [1, 2])
        @test extract_gate_parameters(sCPHASE(1, 2))[1:2] == (:cphase, [π])
        @test extract_gate_parameters(sSWAP(1, 2)) == (:swap, Float64[], [1, 2])
        @test extract_gate_parameters(sXCX(1, 2)) == (:controlled_pauli, Float64[], [1, 2])
        @test extract_gate_parameters(sYCY(1, 2)) == (:controlled_pauli, Float64[], [1, 2])
        @test extract_gate_parameters(sZCZ(1, 2)) == (:controlled_pauli, Float64[], [1, 2])
        @test extract_gate_parameters(sISWAP(1, 2))[1:2] == (:iswap, [π/2])
        @test extract_gate_parameters(sInvISWAP(1, 2))[1:2] == (:inv_iswap, [-π/2])
    end
    
    @testset "Extract Measurement Parameters" begin
        @test extract_gate_parameters(sMX(3, 1))[1] == :measure_x
        @test extract_gate_parameters(sMY(3, 1))[1] == :measure_y
        @test extract_gate_parameters(sMZ(3, 1))[1] == :measure_z
    end
    
    @testset "Extract SparseGate Parameters" begin
        sparse = SparseGate(tCNOT, [2, 5])
        gate_type, params, qubits = extract_gate_parameters(sparse)
        @test qubits == [2, 5]
    end
    
    @testset "Extract PauliMeasurement Parameters" begin
        pauli = P"XYZ"
        pm = PauliMeasurement(pauli, 1)
        gate_type, params, qubits = extract_gate_parameters(pm)
        @test gate_type == :pauli_measurement
        @test qubits == [1, 2, 3]
    end
    
    @testset "T Gate Decomposition" begin
        decomp = get_optimal_gate_decomposition(:T, [π/4], [1])
        @test decomp.gate_type == :T
        @test length(decomp.coefficients) == 2
        @test decomp.target_qubits == [1]
    end
    
    @testset "Phase Gate Decomposition" begin
        decomp = get_optimal_gate_decomposition(:phase, [π/3], [2])
        @test decomp.gate_type == :R_theta
        @test length(decomp.coefficients) == 2
        @test decomp.target_qubits == [2]
    end
    
    @testset "CCZ Gate Decomposition" begin
        decomp = get_optimal_gate_decomposition(:CCZ, Float64[], [1, 2, 3])
        @test decomp.gate_type == :CCZ
        @test length(decomp.coefficients) == 8
        @test decomp.target_qubits == [1, 2, 3]
        
        @test_throws ArgumentError get_optimal_gate_decomposition(:CCZ, Float64[], [1, 2])
    end
    
    @testset "Single Qubit Clifford Decompositions" begin
        for gate_type in [:sqrt_x, :inv_sqrt_x, :sqrt_y, :inv_sqrt_y]
            decomp = get_optimal_gate_decomposition(gate_type, [π/2], [1])
            @test decomp.gate_type == gate_type
            @test length(decomp.coefficients) == 1
            @test decomp.coefficients[1] == ComplexF64(1.0)
        end
        
        for gate_type in [:hadamard, :pauli_x, :pauli_y, :pauli_z, :identity]
            decomp = get_optimal_gate_decomposition(gate_type, Float64[], [1])
            @test decomp.gate_type == gate_type
            @test length(decomp.coefficients) == 1
        end
    end
    
    @testset "Two Qubit Clifford Decompositions" begin
        for gate_type in [:cnot, :cphase, :swap, :controlled_pauli]
            decomp = get_optimal_gate_decomposition(gate_type, Float64[], [1, 2])
            @test decomp.gate_type == gate_type
            @test length(decomp.coefficients) == 1
        end
    end
    
    @testset "Unknown Gate Fallback" begin
        decomp = get_optimal_gate_decomposition(:unknown_two_qubit, Float64[], [1, 2, 3])
        @test decomp.gate_type == :CCZ
        
        decomp = @test_logs (:warn, r"Unknown gate type") get_optimal_gate_decomposition(:totally_unknown, Float64[], [1])
        @test decomp.gate_type == :T
    end
    
    @testset "Pure Clifford Circuit Cost" begin
        circuit = [sHadamard(1), sCNOT(1, 2), sZ(1)]
        cost_info = estimate_simulation_cost(circuit, 0.1)
        
        @test cost_info.estimated_cost == 1
        @test cost_info.total_extent == 1.0
        @test cost_info.non_clifford_gates == 0
        @test isempty(cost_info.individual_extents)
        @test cost_info.scaling_factor == 1.0
    end
    
    @testset "T Gate Cost Estimation" begin
        circuit = [sHadamard(1), TGate(1), TGate(1)]
        cost_info = estimate_simulation_cost(circuit, 0.1)
        
        expected_xi_T = (cos(π/8) + tan(π/8) * sin(π/8))^2
        @test cost_info.non_clifford_gates == 2
        @test cost_info.total_extent ≈ expected_xi_T^2
        @test length(cost_info.individual_extents) == 2
        @test all(xi -> xi ≈ expected_xi_T, cost_info.individual_extents)
    end
    
    @testset "CCZ Gate Cost Estimation" begin
        circuit = [CCZGate(1, 2, 3)]
        cost_info = estimate_simulation_cost(circuit, 0.1)
        
        expected_xi_CCZ = 16.0/9.0
        @test cost_info.non_clifford_gates == 1
        @test cost_info.total_extent ≈ expected_xi_CCZ
        @test cost_info.individual_extents[1] ≈ expected_xi_CCZ
    end
    
    @testset "Mixed Circuit Cost Estimation" begin
        circuit = [
            sHadamard(1), sHadamard(2),
            TGate(1),
            sCNOT(1, 2),
            CCZGate(1, 2, 3),
            sZ(3)
        ]
        cost_info = estimate_simulation_cost(circuit, 0.1)
        
        expected_xi_T = (cos(π/8) + tan(π/8) * sin(π/8))^2
        expected_xi_CCZ = 16.0/9.0
        @test cost_info.non_clifford_gates == 2
        @test cost_info.total_extent ≈ expected_xi_T * expected_xi_CCZ
        @test length(cost_info.individual_extents) == 2
    end
    
    @testset "Cost Scales with Delta" begin
        circuit = [TGate(1)]
        
        cost_small_delta = estimate_simulation_cost(circuit, 0.01)
        cost_large_delta = estimate_simulation_cost(circuit, 0.5)
        
        @test cost_small_delta.estimated_cost > cost_large_delta.estimated_cost
    end
    
    @testset "Valid Construction" begin
        coeffs = [ComplexF64(0.5), ComplexF64(0.3), ComplexF64(0.2)]
        circuits = [
            [sHadamard(1)],
            [sX(1), sZ(1)],
            [sCNOT(1, 2)]
        ]
        
        decomp = CircuitDecomposition(coeffs, circuits, 2)
        
        @test decomp.n_qubits == 2
        @test length(decomp.coefficients) == 3
        @test length(decomp.clifford_circuits) == 3
        @test decomp.l1_norm ≈ sum(abs.(coeffs))
        @test decomp.stabilizer_extent ≈ decomp.l1_norm^2
    end
    
    @testset "Constructor Validation" begin
        coeffs = [ComplexF64(0.5), ComplexF64(0.3)]
        circuits = [[sHadamard(1)]]  
        
        @test_throws DimensionMismatch CircuitDecomposition(coeffs, circuits, 1)
    end
end