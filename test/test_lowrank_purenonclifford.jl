@testitem "PureNonClifford Circuit Simulation" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford
    import QuantumClifford.PureNonClifford:
        get_gate_decomposition,
        CliffordGateDecompositionCache,
        PureGeneralizedStabilizer

    @testset "Pure Clifford Circuit" begin
        circuit_clifford = [sHadamard(1), sCNOT(1,2), sZ(1)]

        result = lrtrajectories(circuit_clifford, 2; trajectories=50, delta=0.1, verbose=false)
        @test result.simulation_cost == 1
        @test result.approximation_error ≈ 0.1
        @test size(lrmeasurements(result), 1) == 50
    end

    @testset "Mixed Circuit" begin
        circuit_mixed = [sHadamard(1), sPhase(1), sCNOT(1,2)]

        result = lrtrajectories(circuit_mixed, 2; trajectories=30, delta=0.2, verbose=false)
        @test result.simulation_cost >= 1
        @test size(lrmeasurements(result), 1) == 30
    end

    @testset "T-Gate Simulation" begin
        circuit_T = [sHadamard(1), sT(1), sHadamard(1)]
        result_T = lrtrajectories(circuit_T, 1; trajectories=100, delta=0.1, verbose=false)

        expected_xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        expected_cost_T = ceil(Int, 1 + expected_xi_T / 0.1^2)

        @test result_T.simulation_cost > 1
        @test result_T.simulation_cost <= expected_cost_T
        @test result_T.approximation_error ≈ 0.1 atol=0.01
    end

    @testset "Multiple T-Gates" begin
        circuit_3T = [sT(1), sT(1), sT(1)]
        result_3T = lrtrajectories(circuit_3T, 1; trajectories=50, delta=0.1, verbose=false)

        expected_xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        expected_xi_3T = expected_xi_T^3
        expected_cost_3T = ceil(Int, 1 + expected_xi_3T / 0.1^2)

        @test result_3T.simulation_cost <= expected_cost_3T
    end

    @testset "CCZ Gate Simulation" begin
        circuit_CCZ = [
            sHadamard(1), sHadamard(2), sHadamard(3),
            sCCZ(1, 2, 3),
            sHadamard(1), sHadamard(2), sHadamard(3)
        ]
        result_CCZ = lrtrajectories(circuit_CCZ, 3; trajectories=50, delta=0.1, verbose=false)

        expected_xi_CCZ = 16.0/9.0
        expected_cost_CCZ = ceil(Int, 1 + expected_xi_CCZ / 0.1^2)

        @test result_CCZ.simulation_cost > 1
        @test result_CCZ.simulation_cost <= expected_cost_CCZ
    end

    @testset "Error Handling" begin
        @test_throws ArgumentError lrtrajectories([], 1)
        @test_throws ArgumentError lrtrajectories([sHadamard(1)], 0)
        @test_throws ArgumentError lrtrajectories([sHadamard(1)], 1; trajectories=-1)
        @test_throws ArgumentError lrtrajectories([sHadamard(1)], 1; delta=1.5)

        circuit_invalid = [sHadamard(5)]
        @test_throws ArgumentError lrtrajectories(circuit_invalid, 2)
    end

    @testset "Cost Estimation" begin
        circuit = [sHadamard(1), sT(1)]
        deltas = [0.5, 0.2, 0.1, 0.05]

        for delta in deltas
            result = lrtrajectories(circuit, 1; trajectories=30, delta=delta, verbose=false)
            xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
            bound = ceil(Int, 1 + xi_T / delta^2)

            @test result.simulation_cost <= bound
            @test result.approximation_error ≈ delta atol=0.02
        end
    end

    @testset "nqubits Interface" begin
        t = sT(1)
        @test nqubits(t) == 1

        ccz = sCCZ(1, 2, 3)
        @test nqubits(ccz) == 3
    end

    @testset "isclifford for non-Clifford gates" begin
        @test isclifford(sT(1)) == false
        @test isclifford(sCCZ(1, 2, 3)) == false
    end

    @testset "AbstractNonCliffordOperator" begin
        @test sT(1) isa QuantumClifford.AbstractNonCliffordOperator
        @test sCCZ(1, 2, 3) isa QuantumClifford.AbstractNonCliffordOperator
        @test sT(1) isa QuantumClifford.AbstractOperation
        @test !(sHadamard(1) isa QuantumClifford.AbstractNonCliffordOperator)
    end

    @testset "stabilizer_extent Trait" begin
        @test stabilizer_extent(sHadamard(1)) == 1.0
        @test stabilizer_extent(sCNOT(1, 2)) == 1.0

        expected_xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        @test stabilizer_extent(sT(1)) ≈ expected_xi_T rtol=0.001

        expected_xi_CCZ = 16.0/9.0
        @test stabilizer_extent(sCCZ(1, 2, 3)) ≈ expected_xi_CCZ rtol=0.001
    end

    @testset "T Gate Decomposition" begin
        decomp = get_gate_decomposition(sT(1))
        @test decomp.gate_type == :T
        @test length(decomp.coefficients) == 2
        @test decomp.target_qubits == [1]
    end

    @testset "CCZ Gate Decomposition" begin
        decomp = get_gate_decomposition(sCCZ(1, 2, 3))
        @test decomp.gate_type == :CCZ
        @test length(decomp.coefficients) == 8
        @test decomp.target_qubits == [1, 2, 3]
    end

    @testset "Pure Clifford Circuit Cost" begin
        circuit = [sHadamard(1), sCNOT(1, 2), sZ(1)]
        cost_info = lrcost(circuit; delta=0.1)

        @test cost_info.estimated_k == 1
        @test cost_info.total_extent == 1.0
        @test cost_info.non_clifford_count == 0
        @test isempty(cost_info.gate_extents)
        @test cost_info.scaling_factor == 1.0
    end

    @testset "T Gate Cost Estimation" begin
        circuit = [sHadamard(1), sT(1), sT(1)]
        cost_info = lrcost(circuit; delta=0.1)

        expected_xi_T = (cos(π/8) + tan(π/8) * sin(π/8))^2
        @test cost_info.non_clifford_count == 2
        @test cost_info.total_extent ≈ expected_xi_T^2
        @test length(cost_info.gate_extents) == 2
        @test all(xi -> xi ≈ expected_xi_T, cost_info.gate_extents)
    end

    @testset "CCZ Gate Cost Estimation" begin
        circuit = [sCCZ(1, 2, 3)]
        cost_info = lrcost(circuit; delta=0.1)

        expected_xi_CCZ = 16.0/9.0
        @test cost_info.non_clifford_count == 1
        @test cost_info.total_extent ≈ expected_xi_CCZ
        @test cost_info.gate_extents[1] ≈ expected_xi_CCZ
    end

    @testset "Mixed Circuit Cost Estimation" begin
        circuit = [
            sHadamard(1), sHadamard(2),
            sT(1),
            sCNOT(1, 2),
            sCCZ(1, 2, 3),
            sZ(3)
        ]
        cost_info = lrcost(circuit; delta=0.1)

        expected_xi_T = (cos(π/8) + tan(π/8) * sin(π/8))^2
        expected_xi_CCZ = 16.0/9.0
        @test cost_info.non_clifford_count == 2
        @test cost_info.total_extent ≈ expected_xi_T * expected_xi_CCZ
        @test length(cost_info.gate_extents) == 2
    end

    @testset "Cost Scales with Delta" begin
        circuit = [sT(1)]

        cost_small_delta = lrcost(circuit; delta=0.01)
        cost_large_delta = lrcost(circuit; delta=0.5)

        @test cost_small_delta.estimated_k > cost_large_delta.estimated_k
    end

    @testset "PureGeneralizedStabilizer construction and apply!" begin
        state = PureGeneralizedStabilizer(2, 0.1)
        @test nqubits(state) == 2
        @test length(state.states) == 1
        @test length(state.coefficients) == 1
        @test state.total_extent == 1.0

        apply!(state, sHadamard(1))
        @test length(state.states) == 1

        apply!(state, sT(1))
        @test length(state.states) == 2
        @test state.total_extent > 1.0
    end
end

@testitem "PureNonClifford Sampling" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford
    import QuantumClifford: AbstractOperation
    import QuantumClifford.PureNonClifford:
        PureGeneralizedStabilizer,
        sample_measurement_outcomes,
        compute_outcome_frequencies

    @testset "Measurement Sampling Tests" begin
        state = PureGeneralizedStabilizer(2, 0.2)
        apply!(state, sHadamard(1))
        apply!(state, sCNOT(1, 2))

        n_samples = 100
        measurements = sample_measurement_outcomes(state, n_samples; verbose=false)

        @test size(measurements, 1) == n_samples
        @test size(measurements, 2) == 2
        @test all(m -> m ∈ [false, true], measurements)

        unique_rows = size(unique(measurements, dims=1), 1)
        @test unique_rows > 0
    end

    @testset "Computational Correctness - Single Qubit" begin
        circuit = [sHadamard(1), sT(1)]

        result = lrtrajectories(circuit, 1; trajectories=1000, delta=0.05, verbose=false)

        measurements = lrmeasurements(result)
        p0 = sum(measurements[:, 1] .== false) / size(measurements, 1)
        p1 = sum(measurements[:, 1] .== true) / size(measurements, 1)

        @test p0 ≈ 0.5 atol=0.05
        @test p1 ≈ 0.5 atol=0.05
        @test result.total_runtime < 30.0
    end

    @testset "Computational Correctness - Entanglement" begin
        circuit = [
            sHadamard(1),
            sCNOT(1,2),
            sT(1)
        ]

        result = lrtrajectories(circuit, 2; trajectories=1000, delta=0.1, verbose=false)

        measurements = lrmeasurements(result)
        freq = compute_outcome_frequencies(measurements)

        p00 = get(freq, BitVector([false, false]), 0.0)
        p11 = get(freq, BitVector([true, true]), 0.0)
        correlation = p00 + p11

        @test correlation > 0.85
        @test result.total_runtime < 60.0
    end

    @testset "Scalability" begin
        qubit_counts = [3, 5, 8]

        for n in qubit_counts
            circuit = AbstractOperation[]
            append!(circuit, [sHadamard(i) for i in 1:n])
            push!(circuit, sT(1))
            append!(circuit, [sHadamard(i) for i in 1:n])

            start = time()
            result = lrtrajectories(circuit, n; trajectories=30, delta=0.2, verbose=false)
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
            sT(1), sT(3), sT(5),
            sCNOT(2,3), sCNOT(4,5),
            sT(2), sT(4)
        ]

        result = lrtrajectories(
            circuit, n_qubits;
            trajectories=200,
            delta=0.1,
            verbose=false
        )

        xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        expected_extent = xi_T^5
        expected_cost_bound = ceil(Int, 1 + expected_extent / 0.1^2)

        @test result.simulation_cost > 1
        @test result.simulation_cost <= expected_cost_bound
        @test result.approximation_error ≈ 0.1 atol=0.02
        @test size(lrmeasurements(result), 1) == 200
        @test result.total_runtime < 120.0
    end
end

@testitem "PureNonClifford Magic States" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford
    using Random
    import QuantumClifford.PureNonClifford:
        get_gate_decomposition,
        decompose_rotation_magic_state,
        MagicStateDecompositionCache,
        sparsify_mixed_destabilizer_decomposition

    @testset "Sparsification (MixedDestabilizer)" begin
        state1 = MixedDestabilizer(S"Z")
        state2 = MixedDestabilizer(S"X")
        state3 = MixedDestabilizer(S"-Z")
        states = [state1, state2, state3]
        coeffs = ComplexF64[0.5, 0.3, 0.24]
        delta = 0.1

        sparse = sparsify_mixed_destabilizer_decomposition(coeffs, states, delta)

        l1_norm = sum(abs.(coeffs))
        theoretical_k_bound = ceil(Int, l1_norm^2 / delta^2)

        @test sparse.k > 0
        @test sparse.k <= theoretical_k_bound
        @test sparse.l1_norm ≈ l1_norm
        @test length(sparse.states) == sparse.k
        @test length(sparse.coefficients) == sparse.k
    end

    @testset "T-Gate Decomposition" begin
        t_decomp = get_gate_decomposition(sT(1))
        expected_xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2

        @test t_decomp.gate_type == :T
        @test length(t_decomp.coefficients) == 2
        @test length(t_decomp.clifford_operations) == 2
        @test t_decomp.stabilizer_extent ≈ expected_xi_T rtol=0.01
    end

    @testset "Rotation Magic State" begin
        theta = π/3
        magic_decomp = decompose_rotation_magic_state(theta)
        expected_xi_R = (cos(theta/2) + tan(π/8)*sin(theta/2))^2

        @test magic_decomp.gate_type == :R_theta
        @test length(magic_decomp.coefficients) == 2
        @test magic_decomp.stabilizer_extent ≈ expected_xi_R rtol=0.01
    end

    @testset "CCZ Gate Decomposition" begin
        ccz_decomp = get_gate_decomposition(sCCZ(1,2,3))
        expected_xi_CCZ = 16.0/9.0

        @test ccz_decomp.gate_type == :CCZ
        @test length(ccz_decomp.coefficients) == 8
        @test length(ccz_decomp.clifford_operations) == 8
        @test ccz_decomp.stabilizer_extent ≈ expected_xi_CCZ rtol=0.01
    end

    @testset "Sparsification with Different Deltas" begin
        states = [MixedDestabilizer(S"X"), MixedDestabilizer(S"Z"),
                  MixedDestabilizer(S"Y"), MixedDestabilizer(S"-X")]
        coeffs = ComplexF64[0.4, 0.3, 0.2, 0.1]

        deltas = [0.5, 0.2, 0.1, 0.05]
        k_values = Int[]

        for delta in deltas
            sparse = sparsify_mixed_destabilizer_decomposition(coeffs, states, delta)
            push!(k_values, sparse.k)
        end

        @test issorted(k_values)
    end

    @testset "Magic State Decomposition Cache Properties" begin
        t_magic = decompose_rotation_magic_state(π/4)
        @test t_magic.stabilizer_fidelity ≈ 1.0 / t_magic.stabilizer_extent rtol=0.01

        ccz_magic = MagicStateDecompositionCache(:CCZ,
            fill(ComplexF64(1/6), 8),
            fill(S"XII IXI IIX", 8))
        @test ccz_magic.stabilizer_extent ≈ 16/9 rtol=0.01
        @test ccz_magic.stabilizer_fidelity ≈ 9/16 rtol=0.01
    end
end

@testitem "PureNonClifford T-Gate Probability Distribution" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford

    @testset "H-T-H Circuit Probability" begin
        circuit = [sHadamard(1), sT(1), sHadamard(1)]

        result = lrtrajectories(circuit, 1; trajectories=10000, delta=0.05, verbose=false)
        measurements = lrmeasurements(result)

        p0 = sum(measurements[:, 1] .== false) / size(measurements, 1)
        expected_p0 = cos(π/8)^2

        @test p0 ≈ expected_p0 atol=0.07

        @test p0 > 0.75
        @test p0 < 0.95
    end

    @testset "T|+⟩ State Distribution" begin
        circuit = [sHadamard(1), sT(1)]

        result = lrtrajectories(circuit, 1; trajectories=5000, delta=0.05, verbose=false)
        measurements = lrmeasurements(result)

        p0 = sum(measurements[:, 1] .== false) / size(measurements, 1)
        p1 = sum(measurements[:, 1] .== true) / size(measurements, 1)

        @test p0 ≈ 0.5 atol=0.08
        @test p1 ≈ 0.5 atol=0.08

        @test p0 > 0.35 && p0 < 0.65
    end

    @testset "Multiple T Gates" begin
        circuit_2T = [sHadamard(1), sT(1), sT(1), sHadamard(1)]

        result = lrtrajectories(circuit_2T, 1; trajectories=5000, delta=0.05, verbose=false)
        measurements = lrmeasurements(result)

        p0 = sum(measurements[:, 1] .== false) / size(measurements, 1)
        expected_p0 = 0.5

        @test p0 ≈ expected_p0 atol=0.08

        @test p0 > 0.35 && p0 < 0.65
    end
end

@testitem "PureNonClifford Bell State with T" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford
    import QuantumClifford.PureNonClifford: compute_outcome_frequencies

    @testset "Bell State Correlations" begin
        circuit = [
            sHadamard(1),
            sCNOT(1, 2),
            sT(1)
        ]

        result = lrtrajectories(circuit, 2; trajectories=5000, delta=0.1, verbose=false)
        measurements = lrmeasurements(result)
        freq = compute_outcome_frequencies(measurements)

        p00 = get(freq, BitVector([false, false]), 0.0)
        p11 = get(freq, BitVector([true, true]), 0.0)
        p01 = get(freq, BitVector([false, true]), 0.0)
        p10 = get(freq, BitVector([true, false]), 0.0)

        correlation = p00 + p11
        anticorrelation = p01 + p10

        @test correlation > 0.9
        @test anticorrelation < 0.1
    end

    @testset "GHZ State with CCZ" begin
        circuit = [
            sHadamard(1), sHadamard(2), sHadamard(3),
            sCCZ(1, 2, 3),
            sHadamard(1), sHadamard(2), sHadamard(3)
        ]

        result = lrtrajectories(circuit, 3; trajectories=2000, delta=0.1, verbose=false)

        @test result.simulation_cost > 1
        @test result.total_extent ≈ 16/9 atol=0.01
        @test size(lrmeasurements(result)) == (2000, 3)
    end
end

@testitem "PureNonClifford Result Display" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford

    @testset "Show Method" begin
        circuit = [sHadamard(1), sT(1)]
        result = lrtrajectories(circuit, 1; trajectories=100, delta=0.1, verbose=false)

        io = IOBuffer()
        show(io, result)
        output = String(take!(io))

        @test occursin("Low-Rank Stabilizer Simulation Results", output)
        @test occursin("Qubits: 1", output)
        @test occursin("Trajectories: 100", output)
        @test occursin("stabilizer terms", output)
    end

    @testset "lrmeasurements Accessor" begin
        circuit = [sHadamard(1), sT(1)]
        result = lrtrajectories(circuit, 1; trajectories=50, delta=0.1, verbose=false)

        m = lrmeasurements(result)

        @test m isa Matrix{Bool}
        @test size(m) == (50, 1)
    end
end

@testitem "PureNonClifford Edge Cases" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford

    @testset "Single Gate Circuits" begin
        result = lrtrajectories([sT(1)], 1; trajectories=50, delta=0.1, verbose=false)
        @test result.simulation_cost > 1

        result = lrtrajectories([sHadamard(1)], 1; trajectories=50, delta=0.1, verbose=false)
        @test result.simulation_cost == 1
    end

    @testset "Large Delta" begin
        circuit = [sHadamard(1), sT(1)]
        result = lrtrajectories(circuit, 1; trajectories=50, delta=0.9, verbose=false)

        @test result.simulation_cost >= 1
        @test result.simulation_cost <= 10
    end

    @testset "Qubit Index Validation" begin
        @test_throws ArgumentError sT(0)
        @test_throws ArgumentError sT(-1)
        @test_throws ArgumentError sCCZ(1, 1, 2)
        @test_throws ArgumentError sCCZ(0, 1, 2)
    end

    @testset "CCZ Gate Constructor" begin
        ccz = sCCZ(3, 1, 2)
        @test ccz.qubits == (1, 2, 3)
    end

    @testset "Auto Qubit Inference" begin
        circuit = [sHadamard(1), sT(3), sCNOT(2, 4)]
        result = lrtrajectories(circuit; trajectories=20, delta=0.2, verbose=false)
        @test result.n_qubits == 4
    end
end

@testitem "PureNonClifford Incremental Sparsification" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford
    using Random

    import QuantumClifford: AbstractOperation
    import QuantumClifford.PureNonClifford:
        sparsify_mixed_destabilizer_decomposition,
        PureGeneralizedStabilizer

    @testset "sparsify_mixed_destabilizer_decomposition" begin
        state1 = MixedDestabilizer(S"Z")
        state2 = MixedDestabilizer(S"X")
        state3 = MixedDestabilizer(S"-Z")
        state4 = MixedDestabilizer(S"-X")

        states = [state1, state2, state3, state4]
        coeffs = ComplexF64[0.5, 0.3, 0.15, 0.05]
        l1_norm = sum(abs, coeffs)

        result_large_delta = sparsify_mixed_destabilizer_decomposition(coeffs, states, 10.0)
        expected_k_large = ceil(Int, l1_norm^2 / 10.0^2)
        @test result_large_delta.k == expected_k_large
        @test result_large_delta.k < length(states)

        result_small_delta = sparsify_mixed_destabilizer_decomposition(coeffs, states, 0.1)
        @test result_small_delta.k == 4
        @test length(result_small_delta.states) == 4
        @test length(result_small_delta.coefficients) == 4

        result_medium_delta = sparsify_mixed_destabilizer_decomposition(coeffs, states, 0.6)
        expected_k_medium = ceil(Int, l1_norm^2 / 0.6^2)
        @test result_medium_delta.k == expected_k_medium
        @test length(result_medium_delta.states) == expected_k_medium
    end

    @testset "Incremental Sparsification Triggers" begin
        n_gates = 20
        circuit = vcat(
            [sHadamard(1)],
            [sT(1) for _ in 1:n_gates],
            [sHadamard(1)]
        )

        delta = 0.7
        result = lrtrajectories(circuit, 1; trajectories=100, delta=delta, verbose=false)

        @test result.simulation_cost < 2^n_gates
        @test result.simulation_cost < 100000

        xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        @test result.total_extent ≈ xi_T^n_gates rtol=0.01
    end

    @testset "Large Circuit Scalability" begin
        n_gates = 24
        circuit = AbstractOperation[sT(1) for _ in 1:n_gates]

        start_time = time()
        result = lrtrajectories(circuit, 1; trajectories=50, delta=0.8, verbose=false)
        elapsed = time() - start_time

        @test elapsed < 120.0
        @test result.simulation_cost < 2^n_gates
        @test result.simulation_cost < 500000

        xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        @test result.total_extent ≈ xi_T^n_gates rtol=0.01
    end

    @testset "Small Circuit No Sparsification" begin
        circuit = AbstractOperation[sT(1), sT(1), sT(1)]

        result = lrtrajectories(circuit, 1; trajectories=50, delta=0.1, verbose=false)

        @test result.simulation_cost == 8
    end

    @testset "Error Budget Distribution" begin
        circuit = AbstractOperation[
            sHadamard(1),
            sT(1), sT(1), sT(1), sT(1),
            sHadamard(1)
        ]

        delta = 0.1
        n_runs = 5
        for _ in 1:n_runs
            result = lrtrajectories(circuit, 1; trajectories=2000, delta=delta, verbose=false)
            measurements = lrmeasurements(result)

            p0 = sum(measurements[:, 1] .== false) / size(measurements, 1)

            @test 0.0 < p0 < 1.0
            @test !isnan(p0)
        end
    end

    @testset "Per-Gate Delta Calculation" begin
        circuit_1T = AbstractOperation[sT(1)]
        circuit_4T = AbstractOperation[sT(1), sT(1), sT(1), sT(1)]

        delta = 0.2

        result_1T = lrtrajectories(circuit_1T, 1; trajectories=50, delta=delta, verbose=false)
        result_4T = lrtrajectories(circuit_4T, 1; trajectories=50, delta=delta, verbose=false)

        @test result_1T.simulation_cost >= 1
        @test result_4T.simulation_cost >= 1

        @test result_1T.approximation_error ≈ delta
        @test result_4T.approximation_error ≈ delta
    end

    @testset "L1 Norm Consistency Check" begin
        circuit = AbstractOperation[sT(1) for _ in 1:10]

        state = lrstate(circuit, 1; delta=0.5)
        l1_final = sum(abs, state.coefficients)

        xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        expected_l1 = sqrt(xi_T^10)

        @test l1_final > 0.5 * expected_l1
        @test l1_final < 2.0 * expected_l1
    end
end

@testitem "PureNonClifford QuantumOptics Comparison" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford: PureGeneralizedStabilizer, lrstate
    using QuantumOpticsBase
    using LinearAlgebra

    b = SpinBasis(1//2)
    T_matrix = Operator(b, ComplexF64[1 0; 0 exp(im*π/4)])
    H_matrix = Operator(b, ComplexF64[1 1; 1 -1]/√2)

    @testset "Single qubit H-T-H circuit" begin
        circuit = [sHadamard(1), sT(1), sHadamard(1)]
        state = lrstate(circuit; delta=0.01)
        ket_qc = Ket(state)

        ket_qo = spinup(b)
        ket_qo = H_matrix * ket_qo
        ket_qo = T_matrix * ket_qo
        ket_qo = H_matrix * ket_qo

        overlap = abs(dagger(ket_qc) * ket_qo)
        @test overlap > 0.85
    end

    @testset "Multiple T gates" begin
        circuit = [sHadamard(1), sT(1), sT(1), sT(1), sHadamard(1)]
        state = lrstate(circuit; delta=0.01)
        ket_qc = Ket(state)

        ket_qo = spinup(b)
        ket_qo = H_matrix * ket_qo
        ket_qo = T_matrix * ket_qo
        ket_qo = T_matrix * ket_qo
        ket_qo = T_matrix * ket_qo
        ket_qo = H_matrix * ket_qo

        overlap = abs(dagger(ket_qc) * ket_qo)
        @test overlap > 0.85
    end
end

@testitem "PureNonClifford Probability Accuracy with Incremental Sparsification" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford

    @testset "H-T^4-H Distribution" begin
        circuit = [sHadamard(1), sT(1), sT(1), sT(1), sT(1), sHadamard(1)]

        result = lrtrajectories(circuit, 1; trajectories=5000, delta=0.05, verbose=false)
        measurements = lrmeasurements(result)

        p1 = sum(measurements[:, 1] .== true) / size(measurements, 1)

        @test p1 > 0.85
    end

    @testset "Multiple T-Gates Accuracy" begin
        circuit = vcat(
            [sHadamard(1)],
            [sT(1) for _ in 1:8],
            [sHadamard(1)]
        )

        result = lrtrajectories(circuit, 1; trajectories=5000, delta=0.1, verbose=false)
        measurements = lrmeasurements(result)

        p0 = sum(measurements[:, 1] .== false) / size(measurements, 1)

        @test p0 > 0.80
    end
end

@testitem "PureNonClifford Incremental Sparsification Correctness" tags=[:non_clifford] begin
    using QuantumClifford
    using QuantumClifford.PureNonClifford
    using Random

    import QuantumClifford: AbstractOperation
    import QuantumClifford.PureNonClifford:
        PureGeneralizedStabilizer,
        compute_outcome_frequencies

    @testset "Probability Accuracy With Sparsification Active" begin
        circuit = AbstractOperation[
            sHadamard(1),
            [sT(1) for _ in 1:16]...,
            sHadamard(1)
        ]

        result = lrtrajectories(circuit, 1; trajectories=3000, delta=0.8, verbose=false)

        @test result.simulation_cost < 65536

        measurements = lrmeasurements(result)
        p0 = sum(measurements[:, 1] .== false) / size(measurements, 1)

        @test p0 > 0.70
    end

    @testset "Multi-Qubit Bell State With Sparsification" begin
        n_T_gates = 16
        circuit = AbstractOperation[
            sHadamard(1),
            sCNOT(1, 2),
            [sT(1) for _ in 1:n_T_gates]...,
        ]

        result = lrtrajectories(circuit, 2; trajectories=2000, delta=0.7, verbose=false)

        @test result.simulation_cost < 2^n_T_gates

        measurements = lrmeasurements(result)
        freq = compute_outcome_frequencies(measurements)

        p00 = get(freq, BitVector([false, false]), 0.0)
        p11 = get(freq, BitVector([true, true]), 0.0)
        p01 = get(freq, BitVector([false, true]), 0.0)
        p10 = get(freq, BitVector([true, false]), 0.0)

        correlation = p00 + p11
        anticorrelation = p01 + p10

        @test correlation > 0.75
        @test anticorrelation < 0.25
    end

    @testset "CCZ With Incremental Sparsification" begin
        circuit = AbstractOperation[
            sHadamard(1), sHadamard(2), sHadamard(3),
            sCCZ(1, 2, 3),
            sCCZ(1, 2, 3),
            sCCZ(1, 2, 3),
            sHadamard(1), sHadamard(2), sHadamard(3)
        ]

        result = lrtrajectories(circuit, 3; trajectories=500, delta=0.9, verbose=false)

        xi_CCZ = 16.0/9.0
        @test result.total_extent ≈ xi_CCZ^3 rtol=0.01

        @test result.simulation_cost >= 1
        @test size(lrmeasurements(result)) == (500, 3)
    end

    @testset "Mixed T and CCZ Gates" begin
        circuit = AbstractOperation[
            sHadamard(1), sHadamard(2), sHadamard(3),
            sT(1),
            sT(2),
            sCCZ(1, 2, 3),
            sT(3),
            sHadamard(1), sHadamard(2), sHadamard(3)
        ]

        result = lrtrajectories(circuit, 3; trajectories=200, delta=0.5, verbose=false)

        xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        xi_CCZ = 16.0/9.0
        expected_extent = xi_T^3 * xi_CCZ
        @test result.total_extent ≈ expected_extent rtol=0.01

        @test result.simulation_cost >= 1
    end

    @testset "Stress Test - Very Large Circuit" begin
        n_gates = 32
        circuit = AbstractOperation[sT(1) for _ in 1:n_gates]

        start_time = time()
        result = lrtrajectories(circuit, 1; trajectories=20, delta=0.95, verbose=false)
        elapsed = time() - start_time

        @test elapsed < 180.0

        @test result.simulation_cost < 1_000_000

        xi_T = (cos(π/8) + tan(π/8)*sin(π/8))^2
        @test result.total_extent ≈ xi_T^n_gates rtol=0.01
    end
end
