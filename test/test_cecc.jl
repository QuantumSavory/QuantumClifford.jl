@testitem "ECC" tags=[:ecc, :ecc_base] begin
    using Nemo
    using Hecke
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractCECC, QECCore

    include("test_cecc_base.jl")

    codes = all_testable_classical_code_instances()

    @testset "correctness checks for classical ECCs" begin
        for code in codes
            H = parity_matrix(code)
            n, k, s  = code_n(code), code_k(code), size(H, 1)
            @test all(col -> any(H[:, col] .!= 0), 1:n)
            @test size(H, 2) == n
            @test size(H, 1) <= size(H, 2)
            H = matrix(GF(2), H)
            computed_rank = rank(H)
            @test n - k == computed_rank
            @test computed_rank <= s && computed_rank <= n
            @test n > 0 && k > 0
            @test k <= n
            rate = k/n
            @test 0 < rate <= 1.0
        end
    end
end

@testitem "BellPair generation" tags=[:ecc, :ecc_base] begin
    using QuantumClifford
    using QuantumClifford.ECC
    using InteractiveUtils: fieldnames

    include("test_cecc_base.jl")

    codes = [BellPairCode(c...) for c in bellpair_circuit_args]

    @test code_n(codes[1]) == 4
    @test code_k(codes[1]) == 2
    @test code_n(codes[2]) == 6
    @test code_k(codes[2]) == 4

    roundtrip = [map(f -> getfield(BellPairCode(c...), f), fieldnames(BellPairCode)) for c in bellpair_circuit_args]
    @test [BellPairCode(c...) for c in roundtrip] == codes
end
