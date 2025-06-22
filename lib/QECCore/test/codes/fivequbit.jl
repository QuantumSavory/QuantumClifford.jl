@testitem "Perfect5" begin
    using Test
    using QECCore

    @testset "Perfect5" begin
        p5 = Perfect5()
        @test parity_matrix(p5) == Bool[1 0 0 1 0 0 1 1 0 0;
                                            0 1 0 0 1 0 0 1 1 0;
                                            1 0 1 0 0 0 0 0 1 1;
                                            0 1 0 1 0 1 0 0 0 1]
        @test code_n(p5) == 5
        @test code_s(p5) == 4
        @test distance(p5) == 3
    end
end