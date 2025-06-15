@testitem "Toric" begin
    using Test
    using QECCore

    @testset "Toric" begin
        c = Toric(2,2)
        pm = Bool[1  0  1  0  1  1  0  0  0  0  0  0  0  0  0  0;
        0  1  0  1  1  1  0  0  0  0  0  0  0  0  0  0;
        1  0  1  0  0  0  1  1  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  1  1  0  0  1  0  1  0;
        0  0  0  0  0  0  0  0  1  1  0  0  0  1  0  1;
        0  0  0  0  0  0  0  0  0  0  1  1  1  0  1  0]
        @test parity_matrix(c) == pm
        @test parity_matrix_x(c) == pm[1:3,1:8]
        @test parity_matrix_z(c) == pm[4:6,9:end]
        @test distance(c) == 2
        @test code_n(c) == 8
        @test code_s(c) == 6
    end
end