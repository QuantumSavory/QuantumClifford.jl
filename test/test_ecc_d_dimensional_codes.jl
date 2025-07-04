@testitem "ECC D-dimensional Surface Code" begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Oscar
        import QECCore: parity_matrix
        using QuantumClifford: stab_looks_good
        using QuantumClifford.ECC: d_dimensional_surface_codes, d_dimensional_toric_codes

        @testset "check stabilizers of D-dimensional surface codes" begin
            @testset "[[L² + (L − 1)², 1, L]] 2D surface code" begin
                for L in 2:5
                    D = 2
                    c = d_dimensional_surface_codes(D, L)
                    code = parity_matrix(c)
                    @test stab_looks_good(code, remove_redundant_rows=true)
                end
            end

            @testset "[L³ + 2L(L − 1)², 1, min(L, L²)]] 3D surface code" begin
                for L in 2:3
                    D = 3
                    c = d_dimensional_surface_codes(D, L)
                    code = parity_matrix(c)
                    @test stab_looks_good(code, remove_redundant_rows=true)
                end
            end

            @testset "[[6L⁴ − 12L³ + 10L² − 4L + 1, 1, L²]] 4D surface code" begin
                L = 2
                D = 4
                c = d_dimensional_surface_codes(D, L)
                code = parity_matrix(c)
                @test stab_looks_good(code, remove_redundant_rows=true)
            end

            @testset "5D surface code" begin
                L = 2
                D = 5
                c = d_dimensional_surface_codes(D, L)
                code = parity_matrix(c)
                @test stab_looks_good(code, remove_redundant_rows=true)
            end
        end

        @testset "check stabilizers of D-dimensional toric codes" begin
            @testset "2D toric code" begin
                for L in 2:5
                    D = 2
                    c = d_dimensional_toric_codes(D, L)
                    code = parity_matrix(c)
                    @test stab_looks_good(code, remove_redundant_rows=true)
                end
            end

            @testset "3D toric code" begin
                for L in 2:3
                    D = 3
                    c = d_dimensional_toric_codes(D, L)
                    code = parity_matrix(c)
                    @test stab_looks_good(code, remove_redundant_rows=true)
                end
            end

            @testset "4D toric code" begin
                L = 2
                D = 4
                c = d_dimensional_toric_codes(D, L)
                code = parity_matrix(c)
                @test stab_looks_good(code, remove_redundant_rows=true)
            end
        end
    end
end
