@testitem "ECC D-dimensional Surface Code" begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Oscar
        using QuantumClifford.ECC: d_dimensional_surface_codes

        @testset "CSS orthogonality check for 2D, 3D, 4D, 5D surface codes" begin
            L = 2
            D = 2
            Hx, Hz′ = d_dimensional_surface_codes(D, L)
            @test iszero(mod.(Hz′*Hx, 2))
            L = 2
            D = 3
            Hx, Hz′ = d_dimensional_surface_codes(D, L)
            @test iszero(mod.(Hz′*Hx, 2))
            L = 2
            D = 4
            Hx, Hz′ = d_dimensional_surface_codes(D, L)
            @test iszero(mod.(Hz′*Hx, 2))
            L = 2
            D = 5
            Hx, Hz′ = d_dimensional_surface_codes(D, L)[2], d_dimensional_surface_codes(D, L)[3]
            @test iszero(mod.(Hz′*Hx, 2))
        end
    end
end
