@testitem "Non-Clifford Circuit Simulation" tags=[:non_clifford] begin
    using QuantumClifford
    import QuantumClifford: AbstractOperation

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
end