@testitem "ECC D-dimensional Surface Code" tags=[:ecc] begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Oscar
        using QECCore
        import HiGHS
        import JuMP
        using Nemo: matrix, GF
        using QuantumClifford: stab_looks_good, stab_to_gf2
        using QuantumClifford.ECC

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
            @testset "2D Toric code" begin
                for L in 2:5
                    D = 2
                    c = DDimensionalToricCode(D, L)
                    code = parity_checks(c)
                    n, k = code_n(code), code_k(code)
                    H = stab_to_gf2(code)
                    mat = matrix(GF(2), H)
                    computed_rank = rank(mat)
                    @test computed_rank == n - k
                    @test code_n(c) == n && code_k(c) == k
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L
                end
            end

            @testset "3D Toric code" begin
                for L in 2:3
                    D = 3
                    c = DDimensionalToricCode(D, L)
                    code = parity_checks(c)
                    n, k = code_n(code), code_k(code)
                    H = stab_to_gf2(code)
                    mat = matrix(GF(2), H)
                    computed_rank = rank(mat)
                    @test computed_rank == n - k
                    @test code_n(c) == n && code_k(c) == k
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == L
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L^2
                end
            end

            @testset "4D Toric code" begin
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
                @test code_n(c) == n && code_k(c) == k
                @test stab_looks_good(code, remove_redundant_rows=true)
                @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
                @test iszero(mod.(metacheck_matrix_x(c)*parity_matrix_x(c), 2))
                @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L^2
                @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L^2
            end
        end
    end
end
