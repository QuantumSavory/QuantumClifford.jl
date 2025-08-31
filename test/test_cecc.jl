@testitem "ECC" tags=[:ecc] begin
    using Nemo
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
