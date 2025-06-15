@testitem "Steane7" begin
    using Test
    using QECCore

    @testset "Steane7" begin
        c = Steane7()
        @test parity_matrix(c) == Bool[0 0 0 1 1 1 1 0 0 0 0 0 0 0;
                                    0 1 1 0 0 1 1 0 0 0 0 0 0 0;
                                    1 0 1 0 1 0 1 0 0 0 0 0 0 0;
                                    0 0 0 0 0 0 0 0 0 0 1 1 1 1;
                                    0 0 0 0 0 0 0 0 1 1 0 0 1 1;
                                    0 0 0 0 0 0 0 1 0 1 0 1 0 1]
        @test code_n(c) == 7
        @test code_s(c) == 6
        @test distance(c) == 3
    end
end