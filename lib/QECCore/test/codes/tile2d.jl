@testitem "Tile 2D" begin
    using Test
    using QECCore
    import QuantumClifford.ECC: parity_checks

    @testset "Tile 2D" begin
        # from table 1 of https://arxiv.org/pdf/2504.09171
        table_I = [
            (288, 8, 3, [(0,0),(2,1),(2,2)], [(0,2),(1,2),(2,0)], 10, 10),  # [[288, 8, 12]]    
            (288, 8, 3, [(0,0),(2,0),(0,1),(0,2)], [(0,0),(0,2),(1,1),(2,2)], 10, 10), # [[288, 8, 14]] 
            (288, 18, 4, [(0,0),(0,3),(2,2),(3,0)], [(0,1),(1,0),(1,1),(3,3)], 9, 9), # [[288, 18, 13]]  
            (512, 18, 4, [(0,0),(0,3),(2,2),(3,0)], [(0,1),(1,0),(1,1),(3,3)], 13, 13)]  # [[512, 18, 19]] 

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
    end
end
