@testitem "ECC D-dimensional Surface Code" tags=[:ecc, :part2] begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Oscar
        using QECCore
        import HiGHS
        import JuMP
        using Nemo: matrix, GF
        using QuantumClifford: stab_looks_good, stab_to_gf2
        using QuantumClifford.ECC
        import Oscar.Generic.MatSpaceElem

        @testset "Check properties of D-dimensional Surface codes" begin
            @testset "[[L² + (L − 1)², 1, L]] 2D surface code" begin
                for L in 2:5
                    D = 2
                    c = DDimensionalSurfaceCode(D, L)
                    code = parity_checks(c)
                    n, k = code_n(code), code_k(code)
                    H = stab_to_gf2(code)
                    mat = matrix(GF(2), H)
                    computed_rank = rank(mat)
                    @test computed_rank == n - k
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test n == L^2 + (L − 1)^2 == code_n(c) && k == 1 == code_k(c)
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L
                end
            end

            @testset "[[L³ + 2L(L − 1)², 1, min(L, L²)]] 3D Surface code" begin
                for L in 2:3
                    D = 3
                    c = DDimensionalSurfaceCode(D, L)
                    code = parity_checks(c)
                    n, k = code_n(code), code_k(code)
                    H = stab_to_gf2(code)
                    mat = matrix(GF(2), H)
                    computed_rank = rank(mat)
                    @test computed_rank == n - k
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test n == L^3 + 2*L*(L − 1)^2 == code_n(c) && k == 1 == code_k(c)
                    @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == L
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L^2
                    dx = distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X))
                    dz = distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z))
                    mindxdz = distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:minXZ))
                    @test min(dx, dz) == mindxdz
                end
            end

            @testset "[[6L⁴ − 12L³ + 10L² − 4L + 1, 1, L²]] 4D Surface code" begin
                # Testing only one instance of 4D codes due to longer execution time.
                L = 2
                D = 4
                c = DDimensionalSurfaceCode(D, L)
                code = parity_checks(c)
                n, k = code_n(code), code_k(code)
                H = stab_to_gf2(code)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k
                @test stab_looks_good(code, remove_redundant_rows=true)
                @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
                @test iszero(mod.(metacheck_matrix_x(c)*parity_matrix_x(c), 2))
                @test n == 6*L^4 − 12*L^3 + 10*L^2 − 4*L + 1 == code_n(c) && k == 1 == code_k(c)
                @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L^2
                @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L^2
            end
        end

        @testset "Check properties of D-dimensional Toric codes" begin
            @testset "[[2L², 2, L]] 2D Toric code" begin
                for L in 2:5
                    D = 2
                    c = DDimensionalToricCode(D, L)
                    code = parity_checks(c)
                    n, k = code_n(code), code_k(code)
                    H = stab_to_gf2(code)
                    mat = matrix(GF(2), H)
                    computed_rank = rank(mat)
                    @test computed_rank == n - k
                    @test code_n(c) == n == 2*L^2 && code_k(c) == k == 2
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L
                end
            end

            @testset "[[3L³, 3, min(L, L²)]] 3D Toric code" begin
                for L in 2:3
                    D = 3
                    c = DDimensionalToricCode(D, L)
                    code = parity_checks(c)
                    n, k = code_n(code), code_k(code)
                    H = stab_to_gf2(code)
                    mat = matrix(GF(2), H)
                    computed_rank = rank(mat)
                    @test computed_rank == n - k
                    @test code_n(c) == n == 3*L^3 && code_k(c) == k == 3
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == L
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L^2
                    dx = distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X))
                    dz = distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z))
                    mindxdz = distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:minXZ))
                    @test min(dx, dz) == mindxdz
                end
            end

            @testset "[[6L⁴, 6, L²]] 4D Toric code" begin
                # Testing only one instance of 4D codes due to longer execution time.
                L = 2
                D = 4
                c = DDimensionalToricCode(D, L)
                code = parity_checks(c)
                n, k = code_n(code), code_k(code)
                H = stab_to_gf2(code)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k
                @test code_n(c) == n == 6*L^4 && code_k(c) == k == 6
                @test stab_looks_good(code, remove_redundant_rows=true)
                @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
                @test iszero(mod.(metacheck_matrix_x(c)*parity_matrix_x(c), 2))
                @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L^2
                @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L^2
            end
        end
    end

    @testset "Homological Product codes of Table III of https://arxiv.org/pdf/2407.18490" begin
        function _circulant_matrix_from_quasi_cyclic_polynomial_matrix(H::MatSpaceElem, l::Int)
            F = GF(2)
            r, n = size(H)
            H_bin = zero_matrix(F, r*l, n*l)
            for i in 1:r, j in 1:n
                cfs = zeros(F, l)
                for (k, c) in enumerate(coefficients(H[i,j]))
                    cfs[k] = c
                end
            circ = hcat([circshift(cfs, s) for s in 0:l-1]...)
            rows = (1:l) .+ (i-1)*l
            cols = (1:l) .+ (j-1)*l
            H_bin[rows, cols] = matrix(F, circ)
            end
            return H_bin
        end
        # [[117, 9, 4]] from Table III
        R, x = polynomial_ring(GF(2), "x")
        l = 3
        H = matrix(R, 2, 3, [x^2 x^2 x^2;
                             x   x^2  0])
        G = matrix(R, 1, 3, [1 x 1+x])
        c = HomologicalProductCode([H,transpose(H)], l)
        code = parity_checks(c)
        n, k = code_n(code), code_k(code)
        H = stab_to_gf2(code)
        mat = matrix(GF(2), H)
        computed_rank = rank(mat)
        @test computed_rank == n - k
        @test n == 117 && k == 9
        @test stab_looks_good(code, remove_redundant_rows=true)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        H = matrix(R, 2, 3, [x^2 x^2 x^2;
                             x   x^2  0])
        G = matrix(R, 1, 3, [1 x 1+x])
        H = _circulant_matrix_from_quasi_cyclic_polynomial_matrix(H, l)
        G = _circulant_matrix_from_quasi_cyclic_polynomial_matrix(G, l)
        @test iszero(H*transpose(G))

        # [[225, 9, 4]] from Table III
        R, x = polynomial_ring(GF(2), "x")
        l = 3
        H = matrix(R, 3, 4, [x^2 x^2 x^2   0;
                             x^2   0 x^2  x^2;
                             x^2 x^2   x  x^2])
        G = matrix(R, 1, 4, [1 (1+x)^2 x^2 (1+x)^2])
        c = HomologicalProductCode([H,transpose(H)], l)
        code = parity_checks(c)
        n, k = code_n(code), code_k(code)
        H = stab_to_gf2(code)
        mat = matrix(GF(2), H)
        computed_rank = rank(mat)
        @test computed_rank == n - k
        @test n == 225 && k == 9
        @test stab_looks_good(code, remove_redundant_rows=true)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6
        H = matrix(R, 3, 4, [x^2 x^2 x^2   0;
                             x^2   0 x^2  x^2;
                             x^2 x^2   x  x^2])
        G = matrix(R, 1, 4, [1 (1+x)^2 x^2 (1+x)^2])
        H = _circulant_matrix_from_quasi_cyclic_polynomial_matrix(H, l)
        G = _circulant_matrix_from_quasi_cyclic_polynomial_matrix(G, l)
        @test iszero(H*transpose(G))

        # [[400, 18, 8]] from Table III
        R, x = polynomial_ring(GF(2), "x")
        l = 4
        H = matrix(R, 3, 4, [x^3 x^3   0 x^3;
                             x^3 x^2 x^3 x^2;
                             x^3 x^3 x^2  0])
        G = matrix(R, 1, 4, [1 1+x+x^2 1+x x+x^2])
        H = _circulant_matrix_from_quasi_cyclic_polynomial_matrix(H, l)
        G = _circulant_matrix_from_quasi_cyclic_polynomial_matrix(G, l)
        @test iszero(H*transpose(G))

        # [[625, 25, 9]]
        R, x = polynomial_ring(GF(2), "x")
        l = 5
        H = matrix(R, 3, 4, [x^4   0 x^4 x^3;
                             0 x^3 x^3 x^4;
                             x^3 x^4   0 x^3])
        G = matrix(R, 1, 4, [1 x+x^2+x^3 1+x^2+x^3 x+x^2])
        H = _circulant_matrix_from_quasi_cyclic_polynomial_matrix(H, l)
        G = _circulant_matrix_from_quasi_cyclic_polynomial_matrix(G, l)
        @test iszero(H*transpose(G))
    end

    @testset "Double Homological Product codes of Table I of https://arxiv.org/pdf/1805.09271" begin
        # [[241, 1, 9]]
        δ = [1 1 0;
             0 1 1]
        c = DoubleHomologicalProductCode(δ)
        code = parity_checks(c)
        n, k = code_n(code), code_k(code)
        H = stab_to_gf2(code)
        mat = matrix(GF(2), H)
        computed_rank = rank(mat)
        @test computed_rank == n - k
        @test n == 241 && k == 1
        @test stab_looks_good(code, remove_redundant_rows=true)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 9
        @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
        @test iszero(mod.(metacheck_matrix_x(c)*parity_matrix_x(c), 2))

        # [[486, 6, 9]]
        δ = [1 1 0;
             0 1 1;
             1 0 1]
        c = DoubleHomologicalProductCode(δ)
        code = parity_checks(c)
        n, k = code_n(code), code_k(code)
        H = stab_to_gf2(code)
        mat = matrix(GF(2), H)
        computed_rank = rank(mat)
        @test computed_rank == n - k
        @test n == 486 && k == 6
        @test stab_looks_good(code, remove_redundant_rows=true)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 9
        @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
        @test iszero(mod.(metacheck_matrix_x(c)*parity_matrix_x(c), 2))

        # [[913, 1, 16]]
        δ = [1 1 0 0;
             0 1 1 0;
             0 0 1 1];
        c = DoubleHomologicalProductCode(δ)
        code = parity_checks(c)
        n, k = code_n(code), code_k(code)
        H = stab_to_gf2(code)
        mat = matrix(GF(2), H)
        computed_rank = rank(mat)
        @test computed_rank == n - k
        @test n == 913 && k == 1
        @test stab_looks_good(code, remove_redundant_rows=true)
        @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
        @test iszero(mod.(metacheck_matrix_x(c)*parity_matrix_x(c), 2))
    end
end
