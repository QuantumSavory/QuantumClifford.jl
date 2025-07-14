@testitem "Surface" begin
    using Test
    using QECCore

    @testset "Surface" begin
        c = Surface(2,2)
        @test parity_matrix(c) == Bool[1 0 1 0 1 0 0 0 0 0;
                                    0 1 0 1 1 0 0 0 0 0;
                                    0 0 0 0 0 1 1 0 0 1;
                                    0 0 0 0 0 0 0 1 1 1]
        @test parity_matrix_x(c) == Bool[1 0 1 0 1;
                                        0 1 0 1 1]
        @test parity_matrix_z(c) == Bool[1 1 0 0 1;
                                        0 0 1 1 1]
        @test distance(c) == 2
        @test code_n(c) == 5
        @test code_s(c) == 4
    end
end
