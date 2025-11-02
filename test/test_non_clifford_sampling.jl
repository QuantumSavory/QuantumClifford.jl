@testitem "Non-Clifford Sampling" tags=[:non_clifford] begin
    using QuantumClifford
    import QuantumClifford: AbstractOperation

    @testset "Measurement Sampling" begin
        state1 = Stabilizer([P"XX"])
        state2 = Stabilizer([P"ZZ"])
        
        states = [state1, state2]
        coeffs = [ComplexF64(0.7), ComplexF64(0.3)]
        
        result = SimulationResult(states, coeffs, 2, 0.1, 1.0)
        
        n_samples = 100
        outcomes = sample_measurement_outcomes(result, n_samples, verbose=false)
        
        @test length(outcomes) == n_samples
        @test all(o -> o.n_qubits == 2, outcomes)
        @test all(o -> all(v -> v ∈ [0,1], o.values), outcomes)
    end

    @testset "Computational Correctness - Single Qubit" begin
        circuit = [sHadamard(1), TGate(1)]
        
        result = simulate_non_clifford_circuit(circuit, 1, samples=1000, delta=0.05, verbose=false)
        
        freq = result.measurement_frequencies
        p0 = get(freq, BitString([0], 1), 0.0)
        p1 = get(freq, BitString([1], 1), 0.0)
        
        @test p0 ≈ 0.5 atol=0.05
        @test p1 ≈ 0.5 atol=0.05
        @test result.total_runtime < 30.0
    end

    @testset "Computational Correctness - Entanglement" begin
        circuit = [
            sHadamard(1),
            sCNOT(1,2),
            TGate(1)
        ]
        
        result = simulate_non_clifford_circuit(circuit, 2, samples=1000, delta=0.1, verbose=false)
        
        freq = result.measurement_frequencies
        p00 = get(freq, BitString([0,0], 2), 0.0)
        p11 = get(freq, BitString([1,1], 2), 0.0)
        correlation = p00 + p11
        
        @test correlation > 0.85
        @test result.total_runtime < 60.0
    end

    @testset "Scalability" begin
        qubit_counts = [3, 5, 8]
        
        for n in qubit_counts
            circuit = AbstractOperation[]
            append!(circuit, [sHadamard(i) for i in 1:n])
            push!(circuit, TGate(1))
            append!(circuit, [sHadamard(i) for i in 1:n])
            
            start = time()
            result = simulate_non_clifford_circuit(circuit, n, samples=30, delta=0.2, verbose=false)
            elapsed = time() - start
            
            @test elapsed < 30.0
            @test result.simulation_cost > 1
        end
    end

    @testset "Complete Pipeline" begin
        n_qubits = 5
        circuit = AbstractOperation[
            [sHadamard(i) for i in 1:n_qubits]...,
            sCNOT(1,2), sCNOT(2,3), sCNOT(3,4), sCNOT(4,5),
            TGate(1), TGate(3), TGate(5),
            sCNOT(2,3), sCNOT(4,5),
            TGate(2), TGate(4)
        ]
        
        result = simulate_non_clifford_circuit(
            circuit, n_qubits, 
            samples=200, 
            delta=0.1, 
            precision=0.01,
            verbose=false
        )
        
        xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        expected_extent = xi_T^5
        expected_cost_bound = ceil(Int, 1 + expected_extent / 0.1^2)
        
        @test result.simulation_cost > 1
        @test result.simulation_cost <= expected_cost_bound
        @test result.approximation_error ≈ 0.1 atol=0.02
        @test length(result.outcomes) == 200
        @test result.total_runtime < 120.0
    end
end