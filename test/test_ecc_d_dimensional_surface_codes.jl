@testitem "ECC D-dimensional Surface Code" begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Oscar
        using QuantumClifford.ECC: d_dimensional_surface_codes

        @testset "CSS orthogonality check D-dimensional surface codes" begin
            @testset "[[L² + (L − 1)², 1, L]] 2D surface code" begin
                for L in 2:5
                    D = 2
                    Hx, Hz′ = d_dimensional_surface_codes(D, L)
                    @test iszero(mod.(Hz′*Hx, 2))
                end
            end

            @testset "[L³ + 2L(L − 1)², 1, min(L, L²)]] 3D surface code" begin
                for L in 2:3
                    D = 3
                    Hx, Hz′ = d_dimensional_surface_codes(D, L)
                    @test iszero(mod.(Hz′*Hx, 2))
                end
            end

            @testset "[[6L⁴ − 12L³ + 10L² − 4L + 1, 1, L²]] 4D surface code" begin
                for L in 2:3
                    D = 4
                    Hx, Hz′ = d_dimensional_surface_codes(D, L)
                    @test iszero(mod.(Hz′*Hx, 2))
                end
            end

            @testset "5D surface code" begin
                for L in 2:3
                    D = 5
                    Hx, Hz′ = d_dimensional_surface_codes(D, L)[2], d_dimensional_surface_codes(D, L)[3]
                    @test iszero(mod.(Hz′*Hx, 2))
                end
            end
        end
    end
end
