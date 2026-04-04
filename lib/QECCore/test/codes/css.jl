@testitem "CSS" begin
    using Test
    using QECCore

    @testset "CSS" begin
        h = QECCore._steane_mat()
        c = CSS(h, h)
        @test parity_matrix(c) == parity_matrix(Steane7())
        @test parity_matrix_x(c) == h
        @test parity_matrix_z(c) == h
    end
end
