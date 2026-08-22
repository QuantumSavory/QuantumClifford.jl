@testitem "Tile2D" begin
    using Test
    using Nemo: matrix, GF, rank
    using QECCore: Tile2D
    using QuantumClifford: stab_looks_good, stab_to_gf2
    using QuantumClifford.ECC: parity_checks, code_n, code_k, parity_matrix_x, parity_matrix_z

    @testset "Constructor" begin
        @test_throws ArgumentError Tile2D(0, [(0,0)], [(0,0)], 1, 1)
        @test_throws ArgumentError Tile2D(2, [(0,0)], [(0,0)], 0, 1)
        @test_throws ArgumentError Tile2D(2, [(0,0)], [(0,0)], 1, 0)
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
        hx = parity_matrix_x(tile)
        hz = parity_matrix_z(tile)
        @test size(hx) == size(hz) == (15, 34)
        @test all(>(0), vec(sum(hx; dims=1)))
        @test all(>(0), vec(sum(hz; dims=1)))
        @test all(>(0), vec(sum(hx; dims=2)))
        @test all(>(0), vec(sum(hz; dims=2)))
        @test all(iseven, hx * transpose(hz))
        @test code_k(tile) == 4

        asymmetric_tile = Tile2D(3, [(0,2), (2,2)], [(2,0), (2,1), (2,2)], 3, 3)
        asymmetric_hx = parity_matrix_x(asymmetric_tile)
        asymmetric_hz = parity_matrix_z(asymmetric_tile)
        @test size(asymmetric_hx) == size(asymmetric_hz) == (15, 34)
        @test all(>(0), vec(sum(asymmetric_hx; dims=1)))
        @test all(>(0), vec(sum(asymmetric_hz; dims=1)))
        @test all(>(0), vec(sum(asymmetric_hx; dims=2)))
        @test all(>(0), vec(sum(asymmetric_hz; dims=2)))
        @test all(iseven, asymmetric_hx * transpose(asymmetric_hz))
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
        nonsquare_hx = parity_matrix_x(nonsquare_tile)
        nonsquare_hz = parity_matrix_z(nonsquare_tile)
        @test nonsquare_hx == expected_x
        @test nonsquare_hz == expected_z
        @test all(iseven, nonsquare_hx * transpose(nonsquare_hz))
        @test code_k(nonsquare_tile) == 2
    end

    @testset "Table I" begin
        # From Table I of https://arxiv.org/pdf/2504.09171
        table_i = [
            (288, 8 , 3, [(0,0), (2,1), (2,2)       ], [(0,2), (1,2), (2,0)       ], 10, 10), # [[288, 8, 12]]
            (288, 8 , 3, [(0,0), (2,0), (0,1), (0,2)], [(0,0), (0,2), (1,1), (2,2)], 10, 10), # [[288, 8, 14]]
            (288, 18, 4, [(0,0), (0,3), (2,2), (3,0)], [(0,1), (1,0), (1,1), (3,3)], 9, 9), # [[288, 18, 13]]
            (512, 18, 4, [(0,0), (0,3), (2,2), (3,0)], [(0,1), (1,0), (1,1), (3,3)], 13, 13)] # [[512, 18, 19]]

        for (n, k, b, horiz, vert, lx, ly) in table_i
            c = Tile2D(b, horiz, vert, lx, ly)
            stab = parity_checks(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            h = stab_to_gf2(stab)
            mat = matrix(GF(2), h)
            computed_rank = rank(mat)
            @test computed_rank == nₛ - kₛ
            @test stab_looks_good(stab, remove_redundant_rows=true)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
        end

        weight_cases = [
            (6, 3, [(0,0), (2,1), (2,2)], [(0,2), (1,2), (2,0)], 10, 10),
            (8, 3, [(0,0), (2,0), (0,1), (0,2)], [(0,0), (0,2), (1,1), (2,2)], 10, 10),
            (8, 4, [(0,0), (0,3), (2,2), (3,0)], [(0,1), (1,0), (1,1), (3,3)], 9, 9),
            (8, 4, [(0,0), (0,3), (2,2), (3,0)], [(0,1), (1,0), (1,1), (3,3)], 13, 13),
        ]
        for (weight, b, horiz_x, vert_x, lx, ly) in weight_cases
            c = Tile2D(b, horiz_x, vert_x, lx, ly)
            @test maximum(sum(parity_matrix_z(c); dims=2)) == weight
            @test maximum(sum(parity_matrix_x(c); dims=2)) == weight
        end
    end

    @testset "Appendix B: Stabilizers of weight-6 confined to 3 by 3 boxes yielding [[288, 8, 12]] codes" begin
        # See https://arxiv.org/pdf/2504.09171

        appendix_b = [
            (288, 8, 3, [(0,0), (0,1), (2,2)], [(0,2), (1,0), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (0,1), (2,2)], [(0,2), (1,2), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (1,0), (2,2)], [(0,1), (0,2), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (1,0), (2,2)], [(0,2), (2,0), (2,1)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (1,2), (2,2)], [(0,1), (0,2), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (1,2), (2,2)], [(0,2), (2,0), (2,1)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (2,1), (2,2)], [(0,2), (1,0), (2,0)], 10, 10),  # [[288, 8, 12]]
            (288, 8, 3, [(0,0), (2,1), (2,2)], [(0,2), (1,2), (2,0)], 10, 10),  # [[288, 8, 12]]
        ]
        for (n, k, b, horiz, vert, lx, ly) in appendix_b
            c = Tile2D(b, horiz, vert, lx, ly)
            stab = parity_checks(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            h = stab_to_gf2(stab)
            mat = matrix(GF(2), h)
            computed_rank = rank(mat)
            @test computed_rank == nₛ - kₛ
            @test stab_looks_good(stab, remove_redundant_rows=true)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
            @test maximum(sum(parity_matrix_z(c); dims=2)) == 6
            @test maximum(sum(parity_matrix_x(c); dims=2)) == 6
        end
    end

    @testset "Appendix C: Stabilizers of weight-8 confined to 3 by 3 boxes yielding [[288, 8, 14]] codes" begin
        # See https://arxiv.org/pdf/2504.09171

        appendix_c = [
            (288, 8, 3, [(0,0), (0,1), (0,2), (2,0)], [(0,0), (0,1), (1,1), (2,2)], 10, 10), # [[288, 8, 14]]
            (288, 8, 3, [(0,0), (0,1), (0,2), (2,0)], [(0,0), (0,2), (1,1), (2,2)], 10, 10), # [[288, 8, 14]]
            (288, 8, 3, [(0,0), (0,1), (0,2), (2,0)], [(0,0), (2,1), (1,2), (2,2)], 10, 10), # [[288, 8, 14]]
            (288, 8, 3, [(0,0), (0,1), (0,2), (2,0)], [(0,1), (1,0), (1,1), (2,2)], 10, 10)  # [[288, 8, 14]]
        ]
        for (n, k, b, horiz, vert, lx, ly) in appendix_c
            c = Tile2D(b, horiz, vert, lx, ly)
            stab = parity_checks(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            h = stab_to_gf2(stab)
            mat = matrix(GF(2), h)
            computed_rank = rank(mat)
            @test computed_rank == nₛ - kₛ
            @test stab_looks_good(stab, remove_redundant_rows=true)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
            @test maximum(sum(parity_matrix_z(c); dims=2)) == 8
            @test maximum(sum(parity_matrix_x(c); dims=2)) == 8
        end
    end
end
