@testitem "ECC D-dimensional Surface Code" tags=[:ecc] begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Oscar
        import QECCore
        import HiGHS
        import JuMP
        using QuantumClifford: stab_looks_good
        using QuantumClifford.ECC

        @testset "check stabilizers of D-dimensional surface codes" begin
            @testset "[[L² + (L − 1)², 1, L]] 2D surface code" begin
                for L in 2:5
                    D = 2
                    c = DDimensionalSurfaceCode(D, L)
                    code = parity_checks(c)
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test code_n(c) == L^2 + (L − 1)^2
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L
                end
            end

            @testset "[[L³ + 2L(L − 1)², 1, min(L, L²)]] 3D surface code" begin
                for L in 2:3
                    D = 3
                    c = DDimensionalSurfaceCode(D, L)
                    code = parity_checks(c)
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test code_n(c) == L^3 + 2*L*(L − 1)^2
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == L
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == L^2
                end
            end

            @testset "[[6L⁴ − 12L³ + 10L² − 4L + 1, 1, L²]] 4D surface code" begin
                L = 2
                D = 4
                c = DDimensionalSurfaceCode(D, L)
                code = parity_checks(c)
                @test stab_looks_good(code, remove_redundant_rows=true)
                @test code_n(c) == 6*L^4 − 12*L^3 + 10*L^2 − 4*L + 1
                @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L^2
            end
        end

        @testset "check stabilizers of D-dimensional toric codes" begin
            @testset "[[L² + (L − 1)², 1, L]] 2D toric code" begin
                for L in 2:5
                    D = 2
                    c = DDimensionalToricCode(D, L)
                    code = parity_checks(c)
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test code_n(c) == L^2 + (L − 1)^2
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L
                end
            end

            @testset "[[L³ + 2L(L − 1)², 1, d]] 3D toric code" begin
                for L in 2:4
                    D = 3
                    c = DDimensionalToricCode(D, L)
                    code = parity_checks(c)
                    @test stab_looks_good(code, remove_redundant_rows=true)
                    @test code_n(c) == L^3 + 2*L*(L − 1)^2
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == L
                end
            end

            @testset "[[6L⁴ − 12L³ + 10L² − 4L + 1, 1, L²]] 4D toric code" begin
                L = 2
                D = 4
                c = DDimensionalToricCode(D, L)
                code = parity_checks(c)
                @test stab_looks_good(code, remove_redundant_rows=true)
                @test code_n(c) == 6*L^4 − 12*L^3 + 10*L^2 − 4*L + 1
                @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == L^2
            end
        end
    end
end
