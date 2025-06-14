@testitem "Shor9" begin
    using Test
    using QECCore

    @testset "Shor9" begin
        c = Shor9()
        @test parity_matrix(c) == Bool[1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
                                    0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0;
                                    0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0;
                                    0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
                                    0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0;
                                    0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0;
                                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0;
                                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1]
        @test code_n(c) == 9
        @test code_s_temp(c) == 8
        @test distance(c) == 3
    end
end