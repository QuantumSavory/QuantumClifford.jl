@testitem "Tile 2D" begin
    using Test
    using Nemo: matrix, GF, rank
    using QECCore: Tile2D
    using QuantumClifford: stab_looks_good, stab_to_gf2
    using QuantumClifford.ECC: parity_checks, code_n, code_k, parity_matrix_x, parity_matrix_z

    @testset "Constructor" begin
        @test_throws ArgumentError Tile2D(0, [(0,0)], [(0,0)], 1, 1)
        @test_throws ArgumentError Tile2D(2, [(-1,0)], [(0,0)], 1, 1)
        @test_throws ArgumentError Tile2D(2, [(0,0)], [(0,2)], 1, 1)

        horiz = [(0,2), (2,2)]
        vert = [(2,0), (2,1), (2,2)]
        tile = Tile2D(3, horiz, vert, 3, 3)
        empty!(horiz)
        empty!(vert)
        @test length(tile.horiz) == 2
        @test length(tile.vert) == 3
    end

    @testset "Parity matrix storage" begin
        tile = Tile2D(2, [(0,0)], [(0,0)], 2, 2)
        @test eltype(parity_matrix_x(tile)) === Bool
        @test eltype(parity_matrix_z(tile)) === Bool
    end

    @testset "Boundary pruning and layout" begin
        tile = Tile2D(3, [(0,2), (1,2), (2,2)], [(2,0), (2,1), (2,2)], 3, 3)
        Hx = parity_matrix_x(tile)
        Hz = parity_matrix_z(tile)
        @test size(Hx) == size(Hz) == (15, 34)
        @test all(>(0), vec(sum(Hx; dims=1)))
        @test all(>(0), vec(sum(Hz; dims=1)))
        @test code_k(tile) == 4

        asymmetric_tile = Tile2D(3, [(0,2), (2,2)], [(2,0), (2,1), (2,2)], 3, 3)
        @test size(parity_matrix_x(asymmetric_tile)) == (15, 34)
        @test size(parity_matrix_z(asymmetric_tile)) == (15, 34)
        @test code_k(asymmetric_tile) == 4

        nonsquare_tile = Tile2D(2, [(0,0), (1,1)], [(0,1), (1,0)], 1, 2)
        expected_x = Bool[
            1 0 0 0 1 0 0 1 0 1 0 0
            0 1 0 0 0 1 0 0 1 0 1 0
            0 0 0 1 0 0 1 0 0 0 0 0
            0 0 1 0 0 0 0 0 0 0 0 1
        ]
        expected_z = Bool[
            0 1 0 1 0 0 1 0 0 0 1 0
            0 0 1 0 1 0 0 1 0 0 0 1
            1 0 0 0 0 0 0 1 0 0 0 0
            0 0 0 0 1 0 0 0 0 1 0 0
            0 1 0 0 0 0 0 0 1 0 0 0
            0 0 0 0 0 1 0 0 0 0 1 0
        ]
        @test parity_matrix_x(nonsquare_tile) == expected_x
        @test parity_matrix_z(nonsquare_tile) == expected_z
        @test code_k(nonsquare_tile) == 2
    end

    @testset "Tile 2D" begin
        # from table 1 of https://arxiv.org/pdf/2504.09171
        table_I = [
            (288, 8 , 3, [(0,0), (2,1), (2,2)       ], [(0,2), (1,2), (2,0)       ], 10, 10), # [[288, 8, 12]]
            (288, 8 , 3, [(0,0), (2,0), (0,1), (0,2)], [(0,0), (0,2), (1,1), (2,2)], 10, 10), # [[288, 8, 14]] 
            (288, 18, 4, [(0,0), (0,3), (2,2), (3,0)], [(0,1), (1,0), (1,1), (3,3)], 9 ,  9), # [[288, 18, 13]]
            (512, 18, 4, [(0,0), (0,3), (2,2), (3,0)], [(0,1), (1,0), (1,1), (3,3)], 13, 13)] # [[512, 18, 19]]

        for (n, k, B, horiz, vert, Lx, Ly) in table_I
            c = Tile2D(B, horiz, vert, Lx, Ly)
            stab = parity_checks(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            H = stab_to_gf2(stab)
            mat = matrix(GF(2), H)
            computed_rank = rank(mat)
            @test computed_rank == nₛ - kₛ
            @test stab_looks_good(stab, remove_redundant_rows=true)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
        end

        # check-weight tests
        # From Table I of https://arxiv.org/pdf/2504.09171v1
        # [[288, 8, 12]] 
        B = 3
        horizX = [(0,0),(2,1),(2,2)]
        vertX = [(0,2),(1,2),(2,0)]
        Lx, Ly = 10, 10
        c = Tile2D(B, horizX, vertX, Lx, Ly)
        @test all(maximum(sum(Matrix(parity_matrix_z(c)), dims=2)) .== 6)
        @test all(maximum(sum(Matrix(parity_matrix_x(c)), dims=2)) .== 6)

        # [[288, 8, 14]]
        B = 3
        horizX = [(0,0),(2,0),(0,1),(0,2)]
        vertX = [(0,0),(0,2),(1,1),(2,2)]
        Lx, Ly = 10, 10
        c = Tile2D(B, horizX, vertX, Lx, Ly)
        @test all(maximum(sum(Matrix(parity_matrix_z(c)), dims=2)) .== 8)
        @test all(maximum(sum(Matrix(parity_matrix_x(c)), dims=2)) .== 8)

        # [[288, 18, 13]]
        B = 4
        horizX = [(0,0),(0,3),(2,2),(3,0)]
        vertX = [(0,1),(1,0),(1,1),(3,3)]
        Lx, Ly = 9, 9
        c = Tile2D(B, horizX, vertX, Lx, Ly)
        @test all(maximum(sum(Matrix(parity_matrix_z(c)), dims=2)) .== 8)
        @test all(maximum(sum(Matrix(parity_matrix_x(c)), dims=2)) .== 8)

        # [[512, 18, 19]]
        B = 4
        horizX = [(0,0),(0,3),(2,2),(3,0)]
        vertX = [(0,1),(1,0),(1,1),(3,3)]
        Lx, Ly = 13, 13
        c = Tile2D(B, horizX, vertX, Lx, Ly)
        @test all(maximum(sum(Matrix(parity_matrix_z(c)), dims=2)) .== 8)
        @test all(maximum(sum(Matrix(parity_matrix_x(c)), dims=2)) .== 8)
    end

    @testset "Appendix B: Stabilizers of weight-6 confined to 3 by 3 boxes yielding [[288, 8, 12]] codes" begin
        # See https://arxiv.org/pdf/2504.09171

        appendix_B = [
            (288, 8, 3, [(0,0), (0,1), (2,2)], [(0,2), (1,0), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (0,1), (2,2)], [(0,2), (1,2), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (1,0), (2,2)], [(0,1), (0,2), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (1,0), (2,2)], [(0,2), (2,0), (2,1)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (1,2), (2,2)], [(0,1), (0,2), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (1,2), (2,2)], [(0,2), (2,0), (2,1)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (2,1), (2,2)], [(0,2), (1,0), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (2,1), (2,2)], [(0,2), (1,2), (2,0)], 10, 10),  # [[288, 8, 12]]
        ]
        for (n, k, B, horiz, vert, Lx, Ly) in appendix_B
            c = Tile2D(B, horiz, vert, Lx, Ly)
            stab = parity_checks(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            H = stab_to_gf2(stab)
            mat = matrix(GF(2), H)
            computed_rank = rank(mat)
            @test computed_rank == nₛ - kₛ
            @test stab_looks_good(stab, remove_redundant_rows=true)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
            @test all(maximum(sum(Matrix(parity_matrix_z(c)), dims=2)) .== 6)
            @test all(maximum(sum(Matrix(parity_matrix_x(c)), dims=2)) .== 6)
        end
    end

    @testset "Appendix C: Stabilizers of weight-8 confined to 3 by 3 boxes yielding [[288, 8, 14]] codes" begin
        # See https://arxiv.org/pdf/2504.09171

        appendix_C = [
            (288, 8, 3, [(0,0), (0,1), (0,2), (2,0)], [(0,0), (0,1), (1,1), (2,2)], 10, 10), # [[288, 8, 14]]
            (288, 8, 3, [(0,0), (0,1), (0,2), (2,0)], [(0,0), (0,2), (1,1), (2,2)], 10, 10), # [[288, 8, 14]]
            (288, 8, 3, [(0,0), (0,1), (0,2), (2,0)], [(0,0), (2,1), (1,2), (2,2)], 10, 10), # [[288, 8, 14]]
            (288, 8, 3, [(0,0), (0,1), (0,2), (2,0)], [(0,1), (1,0), (1,1), (2,2)], 10, 10)  # [[288, 8, 14]]
        ]
        for (n, k, B, horiz, vert, Lx, Ly) in appendix_C
            c = Tile2D(B, horiz, vert, Lx, Ly)
            stab = parity_checks(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            H = stab_to_gf2(stab)
            mat = matrix(GF(2), H)
            computed_rank = rank(mat)
            @test computed_rank == nₛ - kₛ
            @test stab_looks_good(stab, remove_redundant_rows=true)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
            @test all(maximum(sum(Matrix(parity_matrix_z(c)), dims=2)) .== 8)
            @test all(maximum(sum(Matrix(parity_matrix_x(c)), dims=2)) .== 8)
        end
    end
end
