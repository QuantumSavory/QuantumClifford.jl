@testitem "ECC Mirror Code" tags=[:ecc, :ecc_bespoke_checks, :oscar_required] begin
    using Oscar
    using QECCore
    import HiGHS
    import JuMP
    using Nemo: matrix, GF
    using QuantumClifford: stab_looks_good, stab_to_gf2
    using QuantumClifford.ECC

    @testset "Mirror Codes" begin

        table_i_product_of_three_groups = [
            (30 , 4 , 5 , 2, 3, [2, 3,  5], [(0, 0, 0), (0, 0, 1)], [(1, 0, 0), (1, 1, 0), (1, 2, 2)]),
            (36 , 4 , 6 , 2, 3, [4, 3,  3], [(0, 0, 0), (2, 0, 1)], [(1, 0, 0), (1, 1, 0), (1, 2, 1)]),
            (42 , 6 , 6 , 2, 4, [2, 3,  7], [(0, 0, 0), (0, 0, 1)], [(1, 0, 0), (1, 0, 2), (1, 1, 0), (1, 1, 3)]),
            (70 , 10, 6 , 2, 4, [2, 5,  7], [(0, 0, 0), (0, 0, 1)], [(1, 0, 0), (1, 0, 2), (1, 1, 0), (1, 1, 3)]),
            (80 , 6 , 10, 2, 4, [2, 8,  5], [(0, 0, 0), (1, 0, 1)], [(0, 1, 0), (0, 3, 0), (1, 5, 1), (1, 7, 2)]),
            (90 , 10, 7 , 2, 4, [2, 9,  5], [(0, 0, 0), (0, 1, 0)], [(1, 0, 0), (1, 0, 1), (1, 2, 0), (1, 4, 1)]),
            (126, 14, 7 , 2, 4, [2, 9,  7], [(0, 0, 0), (0, 1, 0)], [(1, 0, 0), (1, 0, 1), (1, 2, 0), (1, 4, 1)]),
            (132, 8 , 11, 2, 4, [4, 3, 11], [(0, 0, 0), (0, 0, 1)], [(1, 0, 0), (1, 1, 1), (3, 0, 2), (3, 1, 5)]),
            (174, 6 , 16, 2, 4, [3, 3, 29], [(0, 0, 0), (0, 0, 1)], [(1, 0, 0), (1, 0, 4), (1, 1, 7), (1, 1, 18)]),
        ]

        for (n, k, d, wx, wz, orders, A, B) in table_i_product_of_three_groups
            G = abelian_group(orders)
            c = Mirror(G, A, B, true)
            stab = parity_checks(c)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(stab) - code_k(stab)
            @test n == code_n(stab)
            @test k == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
        end
    end
end
