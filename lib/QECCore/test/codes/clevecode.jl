@testitem "Cleve8" begin
    using Test
    using QECCore

    @testset "Cleve8" begin
        c = Cleve8()
        @test parity_matrix(c) == Bool[1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
                                    0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
                                    1 0 1 0 0 1 0 1 0 0 0 0 1 1 1 1;
                                    1 0 1 0 1 0 1 0 0 0 1 1 0 0 1 1;
                                    1 0 0 1 0 1 1 0 0 1 0 1 0 1 0 1]
        @test code_n(c) == 8
        @test code_s_temp(c) == 5
    end
end