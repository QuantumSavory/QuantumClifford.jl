@testitem "ECC Lacross" begin
    using JuMP
    using HiGHS
    using Hecke
    using QuantumClifford
    using QuantumClifford: stab_looks_good, gf2_row_echelon_with_pivots!
    using QuantumClifford.ECC
    using Nemo: matrix, GF
    using QECCore
    
    function _scratch_matrix(n, k, coeffs)
        first_row = zeros(Int, n)
        first_row[1] = 1 # constant term (x⁰)
        for i in 1:k
            first_row[i+1] = coeffs[i] # coefficients for x¹, x², ..., xᵏ
        end
        H = zeros(Int, n, n)
        H[1, :] = first_row
        for i in 2:n
            H[i, :] = circshift(H[i-1, :], 1)
        end
        return H
    end

    function _hecke_circulant_matrix(n, h)
        R = parent(h)
        x = gen(R)
        _, proj = residue_ring(R, R(x)^n-1)
        h = proj(h)
        lifted_h = lift(h)
        coeffs = Int[lift(ZZ, coeff(lifted_h, i)) for i in 0:n-1]
        H = zero_matrix(GF(2), n, n)
        for i in 1:n
            for j in 1:n
                H[i, j] = coeffs[mod1(j-i+1, n)]
            end
        end
        H = [Int(lift(ZZ, H[i,j])) for i in 1:nrows(H), j in 1:ncols(H)]
        return H
    end

    @testset "circulant matrix consistency check" begin
        for _ in 1:20
            n = rand(1:20)
            k = rand(0:n-1)
            coeffs = rand([0, 1], k)
            @testset "(n=$n, k=$k)" begin
                H = _scratch_matrix(n, k, coeffs)
                R, x = polynomial_ring(GF(2), "x")
                h = isempty(coeffs) ? x^0 : sum(coeffs[i]*x^i for i in 1:k) + x^0
                Hₕ = _hecke_circulant_matrix(n, h)
                @test H == Hₕ
            end
        end
    end

    @testset "Reproduce Figure 3 of [pecorari2025high](@cite)" begin
        @testset "Reproduce Figure 3(a) of [pecorari2025high](@cite)" begin
            # [[52, 4, 4]]
            n = 6
            F = GF(2)
            R, x = polynomial_ring(F, "x")
            h = 1 + x + x^2;
            full_rank = true
            c = Lacross(n, h, full_rank)
            mat = matrix(GF(2), stab_to_gf2(parity_checks(c)))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == 52 && code_k(c) == 4
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
            @test stab_looks_good(parity_checks(c), remove_redundant_rows=true) == true

            # [[130, 4, 6]]
            n = 9
            full_rank = true
            c = Lacross(n, h, full_rank)
            mat = matrix(GF(2), stab_to_gf2(parity_checks(c)))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == 130 && code_k(c) == 4
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6
            @test stab_looks_good(parity_checks(c), remove_redundant_rows=true) == true

            # [[244, 4, 8]]
            n = 12
            full_rank = true
            c = Lacross(n, h, full_rank)
            mat = matrix(GF(2), stab_to_gf2(parity_checks(c)))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == 244 && code_k(c) == 4
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8
            @test stab_looks_good(parity_checks(c), remove_redundant_rows=true) == true
        end

        @testset "Reproduce Figure 3(b) of [pecorari2025high](@cite)" begin
            # [[65, 9, 4]]
            n = 7
            F = GF(2)
            R, x = polynomial_ring(F, "x")
            h = 1 + x + x^3;
            full_rank = true
            c = Lacross(n, h, full_rank)
            mat = matrix(GF(2), stab_to_gf2(parity_checks(c)))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == 65 && code_k(c) == 9
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
            @test stab_looks_good(parity_checks(c), remove_redundant_rows=true) == true

            # [[98, 18, 4]]
            full_rank = false
            c = Lacross(n, h, full_rank)
            mat = matrix(GF(2), stab_to_gf2(parity_checks(c)))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == 98 && code_k(c) == 18
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
            @test stab_looks_good(parity_checks(c), remove_redundant_rows=true) == true
        end

        @testset "cross-checks from pg.4 of https://arxiv.org/pdf/2404.13010" begin
            # [[117, 9, d]]
            n = 9
            F = GF(2)
            R, x = polynomial_ring(F, "x")
            h = 1 + x^3;
            full_rank = true # corresponds to code with open boundary conditions
            c = Lacross(n, h, full_rank)
            mat = matrix(GF(2), stab_to_gf2(parity_checks(c)))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == 117 && code_k(c) == 9
            @test stab_looks_good(parity_checks(c), remove_redundant_rows=true) == true

            # [[162, 18, d]]
            full_rank = false # corresponds to code with periodic boundary conditions
            c = Lacross(n, h, full_rank)
            mat = matrix(GF(2), stab_to_gf2(parity_checks(c)))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == 162 && code_k(c) == 18
            @test stab_looks_good(parity_checks(c), remove_redundant_rows=true) == true
        end
    end
end
