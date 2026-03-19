@testitem "ECC ZSZ" tags=[:ecc, :ecc_bespoke_checks, :oscar_required] begin
    using Oscar
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: DistanceMIPAlgorithm, code_n, code_k, distance, ZSZ, parity_matrix_x, parity_matrix_z, parity_checks
    using QuantumClifford: stab_looks_good, stab_to_gf2

    # all 13 codes from table 1 of Guo et al. (arXiv:2507.21396)
    # format: (name, l, m, q, A, B, expected_k)
    all_table1_codes = [
        ("ZSZ80",    5,  8,  2,  [(0,0), (4,4), (4,1)],   [(0,0), (3,0), (2,7)],   2),
        ("ZSZ108",   3,  18, 2,  [(0,0), (1,0), (0,3)],   [(0,0), (1,1), (2,12)],  2),
        ("ZSZ144-3", 12, 6,  5,  [(0,1), (2,2), (8,3)],   [(4,0), (6,4), (1,5)],   12),
        ("ZSZ160",   5,  16, 2,  [(0,0), (0,9), (4,14)],  [(0,0), (1,13), (2,13)], 2),
        ("ZSZ162",   27, 3,  10, [(8,0), (18,1), (25,2)], [(21,0), (14,1), (16,2)],8),
        ("ZSZ180",   3,  30, 2,  [(0,0), (2,0), (0,3)],   [(0,0), (1,18), (1,19)], 2),
        ("ZSZ288-1", 24, 6,  5,  [(1,1), (3,2), (5,3)],   [(10,3), (23,4), (21,4)],12),
        ("ZSZ288-2", 24, 6,  5,  [(0,0), (14,0), (3,3)],  [(0,0), (0,1), (8,5)],   12),
        ("ZSZ360-1", 30, 6,  19, [(0,0), (8,2), (26,4)],  [(0,0), (21,4), (5,5)],  16),
        ("ZSZ360-2", 30, 6,  19, [(0,0), (9,1), (26,4)],  [(0,0), (4,5), (16,3)],  20),
        ("ZSZ360-3", 30, 6,  11, [(12,1), (21,4), (4,5)], [(4,3), (20,3), (22,4)], 12),
        ("ZSZ540",   45, 6,  19, [(0,0), (34,0), (4,1)],  [(0,0), (38,0), (2,5)],  16),
        ("ZSZ756",   42, 9,  25, [(0,0), (4,4), (6,3)],   [(0,0), (36,5), (21,7)], 24)
    ]

    @testset "Construction and code_n for all Table 1 codes" begin
        for (name, l, m, q, A, B, expected_k) in all_table1_codes
            c = ZSZ(l, m, q, A, B)
            @test code_n(c) == 2 * l * m
        end
    end

    @testset "Full parity checks for all Table 1 codes" begin
        for (name, l, m, q, A, B, expected_k) in all_table1_codes
            c = ZSZ(l, m, q, A, B)
            expected_n = 2 * l * m

            stab = parity_checks(c)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == expected_n == code_n(stab)
            @test code_k(c) == expected_k == code_k(stab)

            Hx = parity_matrix_x(c)
            Hz = parity_matrix_z(c)

            @test all(sum(Hx, dims=1) .== 3)
            @test all(sum(Hx, dims=2) .== 6)
            @test all(sum(Hz, dims=1) .== 3)
            @test all(sum(Hz, dims=2) .== 6)
        end
    end

    @testset "Invalid ZSZ Parameters" begin
        # q^m mod l = 2^3 mod 5 = 3 ≠ 1, so this should be rejected
        @test_throws ArgumentError ZSZ(5, 3, 2, [(0,0)], [(0,0)])
    end
end