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
            n, k, s = code_n(code), code_k(code), size(H, 1)
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

    include("test_ecc_base.jl")

    # constructor / dimension tests
    codes = [BellPairCode(c...) for c in bellpair_circuit_args]

    @test code_n(codes[1]) == 4
    @test code_k(codes[1]) == 2
    @test code_n(codes[2]) == 6
    @test code_k(codes[2]) == 4

    # round-trip equality test (requires Base.== on BellPairCode)
    roundtrip = [map(f -> getfield(BellPairCode(c...), f), fieldnames(BellPairCode)) for c in bellpair_circuit_args]
    @test [BellPairCode(c...) for c in roundtrip] == codes

    # test parity_checks with a real non-trivial circuit, not just empty ones
    rng = MersenneTwister(42)
    for (n, k) in [(2, 1), (3, 1), (4, 2)]
        ngates = 2 * n
        code = random_all_to_all_bellpair_code(rng, n, ngates, k)
        P = parity_checks(code)
        @test size(P, 2) == code_n(code)           # correct number of physical qubits
        @test size(P, 1) == code_n(code) - code_k(code)  # n - k stabilisers
    end

    # brickwork variant with real circuit
    for (n, k) in [(2, 1), (4, 2)]
        nlayers = 3
        code = random_brickwork_bellpair_code(rng, n, nlayers, k)
        P = parity_checks(code)
        @test size(P, 2) == code_n(code)
        @test size(P, 1) == code_n(code) - code_k(code)
    end

    # encode_pairs ordering: same slots in different order should produce equal codes
    c1 = BellPairCode(3, QuantumClifford.AbstractOperation[], [1, 3])
    c2 = BellPairCode(3, QuantumClifford.AbstractOperation[], [3, 1])
    @test c1 == c2
end
