@testitem "Non-Clifford Magic States" tags=[:non_clifford] begin
    using QuantumClifford
    using Random

    @testset "Sparsification Lemma" begin
        states = [S"X", S"Z", S"Y"]
        coeffs = [ComplexF64(0.5), ComplexF64(0.3), ComplexF64(0.24)]
        delta = 0.1
        
        sparse = sparsify_stabilizer_decomposition(coeffs, states, delta)
        
        l1_norm = sum(abs.(coeffs))
        theoretical_k_bound = ceil(Int, 1 + l1_norm^2 / delta^2)
        
        @test sparse.k > 0
        @test sparse.k <= theoretical_k_bound
        @test sparse.approximation_error_bound == delta
        @test sparse.original_l1_norm ≈ l1_norm
        @test length(sparse.states) == sparse.k
        @test length(sparse.coefficients) == sparse.k
    end

    @testset "T-Gate Decomposition" begin
        t_decomp = decompose_T_gate(1)
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
        ccz_decomp = decompose_CCZ_gate([1,2,3])
        expected_xi_CCZ = 16.0/9.0
        
        @test ccz_decomp.gate_type == :CCZ
        @test length(ccz_decomp.coefficients) == 8
        @test length(ccz_decomp.clifford_operations) == 8
        @test ccz_decomp.stabilizer_extent ≈ expected_xi_CCZ rtol=0.01
    end

    @testset "Norm Estimation Tests" begin
        states = [S"X", S"Z", S"Y"]
        coeffs = [ComplexF64(0.5), ComplexF64(0.3), ComplexF64(0.24)]
        sparse = sparsify_stabilizer_decomposition(coeffs, states, 0.1)
        
        estimates = Float64[]
        for _ in 1:3
            est = estimate_sparsification_quality(sparse)
            push!(estimates, est.expected_norm)
        end
        
        @test all(e -> e > 0, estimates)
        @test length(estimates) == 3
    end
end