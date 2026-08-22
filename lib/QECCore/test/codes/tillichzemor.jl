@testitem "Quantum Tillich-Zemor materialized ownership" tags=[:qec_determinism] begin
    using Test
    using SparseArrays
    using QECCore

    deterministic = TillichZemor(4, 3, 3)
    source_Hx, source_Hz = parity_matrix_xz(deterministic)
    expected_Hx, expected_Hz = copy(source_Hx), copy(source_Hz)
    materialized = TillichZemor(4, 3, 3, (source_Hx, source_Hz))

    source_Hx .= 0
    source_Hz .= 0
    @test parity_matrix_xz(materialized) == (expected_Hx, expected_Hz)
    @test parity_matrix_x(materialized) isa Matrix{Int}
    @test parity_matrix_z(materialized) isa Matrix{Int}

    returned_Hx, returned_Hz = parity_matrix_xz(materialized)
    returned_Hx .= 0
    returned_Hz .= 0
    @test parity_matrix_xz(materialized) == (expected_Hx, expected_Hz)
    @test all(iszero, mod.(expected_Hx * transpose(expected_Hz), 2))

    sparse_Hx, sparse_Hz = sparse(expected_Hx), sparse(expected_Hz)
    sparse_code = TillichZemor(4, 3, 3, (sparse_Hx, sparse_Hz))
    sparse_Hx .= 0
    sparse_Hz .= 0
    @test parity_matrix_x(sparse_code) isa SparseMatrixCSC{Int}
    @test parity_matrix_z(sparse_code) isa SparseMatrixCSC{Int}
    @test parity_matrix_xz(sparse_code) == (sparse(expected_Hx), sparse(expected_Hz))

    reconstructed = TillichZemor(4, 3, 3, (copy(expected_Hx), copy(expected_Hz)))
    @test parity_matrix(reconstructed) == parity_matrix(materialized)
end

@testitem "Quantum Tillich-Zemor" begin

    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using Nemo: matrix, GF
    using QECCore.LinearAlgebra
    using QECCore
    using QECCore: _create_circulant_matrix, _create_matrix_M_deterministic, random_TillichZemor_code
    using Random
    using Random: AbstractRNG

    function _theoretical_code_k(c::TillichZemor)
        C = _create_circulant_matrix(c.m)
        M = _create_matrix_M_deterministic(c.m, c.n, c.r)
        H = hcat(C, M)
        H_gf2 = matrix(GF(2), H)
        Ht_gf2 = transpose(H_gf2)
        k = c.n - rank(H_gf2)
        kT = c.m - rank(Ht_gf2) # c.m == size(H, 1)
        k_q = k^2 + kT^2
        return k_q
    end

    @testset "Testing Quantum Tillich-Zemor properties" begin
        for n in 4:2:20
            for m in 3:10
                for r in 3:5
                    if m ≥ r && (n - m)*r ≥ m
                        c = TillichZemor(n, m, r)
                        stab = parity_checks(c)
                        nₛ, kₛ = code_n(stab), code_k(stab)
                        H = stab_to_gf2(stab)
                        mat = matrix(GF(2), H)
                        computed_rank = rank(mat)
                        @test computed_rank == nₛ - kₛ
                        @test stab_looks_good(stab, remove_redundant_rows=true)
                        # test the derivation of [[N, K, _ ]] formules for these novel codes.
                        @test code_n(c) == nₛ && code_k(c) == kₛ == _theoretical_code_k(c)
                    end
                end
            end
        end
    end

    rng = MersenneTwister(1234)
    @testset "Testing random Quantum Tillich-Zemor properties" begin
        for n in 4:2:20
            for m in 3:10
                for r in 3:5
                    if m ≥ r && (n - m)*r ≥ m
                        c = random_TillichZemor_code(rng, n, m, r)
                        stab = parity_checks(c)
                        nₛ, kₛ = code_n(stab), code_k(stab)
                        H = stab_to_gf2(stab)
                        mat = matrix(GF(2), H)
                        computed_rank = rank(mat)
                        @test computed_rank == nₛ - kₛ
                        @test stab_looks_good(stab, remove_redundant_rows=true)
                        @test code_n(c) == nₛ && code_k(c) == kₛ
                    end
                end
            end
        end
    end
end
