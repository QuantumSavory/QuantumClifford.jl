@testitem "ECC GHP Code" tags=[:ecc] begin
    using Hecke
    using QuantumClifford
    using QuantumClifford: stab_looks_good, gf2_row_echelon_with_pivots!
    using QuantumClifford.ECC
    using Nemo: matrix, GF
    using QECCore
    using QuantumClifford.ECC: random_qc_ghp_code_matrix_A

    @testset "Appendix B of [panteleev2021degenerate](@cite)" begin
        # [[882, 24, 18 ≤ d ≤ 24]]
        F = GF(2)
        R, x = polynomial_ring(F, "x")
        n = 7
        l = 63
        S, _ =  quo(R, x^l - 1)
        A = matrix(S, n, n,
            [x^27  0     0     0     0     1     x^54
             x^54  x^27  0     0     0     0     1
             1     x^54  x^27  0     0     0     0
             0     1     x^54  x^27  0     0     0
             0     0     1     x^54  x^27  0     0
             0     0     0     1     x^54  x^27  0
             0     0     0     0     1     x^54  x^27])
        b = S(1 + x + x^6)
        c = GeneralizedHyperGraphProductCode(A, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 882 == code_n(stab) && code_k(c) == 24 == code_k(stab)
        @test stab_looks_good(stab, remove_redundant_rows=true)
        @test all(sum(parity_matrix_x(c), dims=1) .== 3)
        @test all(sum(parity_matrix_x(c), dims=2) .== 6)
        @test all(sum(parity_matrix_z(c), dims=1) .== 3)
        @test all(sum(parity_matrix_z(c), dims=2) .== 6)

        # [[882, 48, 16]]
        F = GF(2)
        R, x = polynomial_ring(F, "x")
        n = 7
        l = 63
        S, _ =  quo(R, x^l - 1)
        A = matrix(S, n, n,
            [x^27   0     0     1     x^18  x^27  1
             1      x^27  0     0     1     x^18  x^27
             x^27   1     x^27  0     0     1     x^18
             x^18   x^27  1     x^27  0     0     1
             1      x^18  x^27  1     x^27  0     0
             0      1     x^18  x^27  1     x^27  0
             0      0     1     x^18  x^27  1     x^27])
        b = S(1 + x + x^6)
        c = GeneralizedHyperGraphProductCode(A, b, l)
        stab = parity_checks(c);
        mat = matrix(GF(2), stab_to_gf2(stab));
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 882 == code_n(stab) && code_k(c) == 48 == code_k(stab)
        @test stab_looks_good(stab, remove_redundant_rows=true)
        @test all(sum(parity_matrix_x(c), dims=2) .== 8)
        @test all(sum(parity_matrix_z(c), dims=2) .== 8)
        @test count(==(3), vec(sum(parity_matrix_x(c), dims=1))) == size(parity_matrix_x(c), 2) ÷ 2
        @test count(==(5), vec(sum(parity_matrix_x(c), dims=1))) == size(parity_matrix_x(c), 2) ÷ 2
        @test count(==(3), vec(sum(parity_matrix_z(c), dims=1))) == size(parity_matrix_z(c), 2) ÷ 2
        @test count(==(5), vec(sum(parity_matrix_z(c), dims=1))) == size(parity_matrix_z(c), 2) ÷ 2

        # [[1270, 28, 16≤ d ≤ 46]]
        n = 5
        l = 127
        S, _ =  quo(R, x^l - 1)
        A = matrix(S, n, n,
            [1     0     x^51  x^52  0
             0     1     0     x^111 x^20
             1     0     x^98  0     x^122
             1     x^80  0     x^119 0 
             0     1     x^5   0     x^106])
        b = S(1 + x + x^7)
        c = GeneralizedHyperGraphProductCode(A, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 1270 == code_n(stab) && code_k(c) == 28 == code_k(stab)
        @test stab_looks_good(stab, remove_redundant_rows=true)
        @test all(sum(parity_matrix_x(c), dims=1) .== 3)
        @test all(sum(parity_matrix_x(c), dims=2) .== 6)
        @test all(sum(parity_matrix_z(c), dims=1) .== 3)
        @test all(sum(parity_matrix_z(c), dims=2) .== 6)
    end
    
    @testset "random quasi-cyclic A matrix tests for quasi-cyclic GHP codes" begin 
        F = GF(2)
        R, x = polynomial_ring(F, "x")
        l = 63
        n = 7
        w = 3
        S, _ =  quo(R, x^l - 1)
        b = S(1 + x + x^6)
        for _ in 1:5
            for n in 5:7
                b_w = degree(lift(b))
                A = random_qc_ghp_code_matrix_A(S, b,  n, w, l)
                c = GeneralizedHyperGraphProductCode(A, b, l)
                stab = parity_checks(c)
                mat = matrix(GF(2), stab_to_gf2(stab))
                computed_rank = rank(mat)
                @test computed_rank == code_n(c) - code_k(c)
                @test stab_looks_good(stab, remove_redundant_rows=true)
                @test all(sum(parity_matrix_x(c), dims=1) .== w)
                @test all(sum(parity_matrix_z(c), dims=2) .== b_w)
                @test all(sum(parity_matrix_x(c), dims=1) .== w)
                @test all(sum(parity_matrix_z(c), dims=2) .== b_w)
            end
        end
    end
end
