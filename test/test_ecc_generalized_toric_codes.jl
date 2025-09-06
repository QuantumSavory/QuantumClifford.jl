@testitem "ECC Generalized Toric Code" tags=[:ecc] begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Oscar
        using QECCore
        import HiGHS
        import JuMP
        using Nemo: matrix, GF
        using QuantumClifford: stab_looks_good, stab_to_gf2
        using QuantumClifford.ECC

        @testset "Generalized Bicycle Codes from Tables 1-IV of https://arxiv.org/abs/2503.03827" begin
            R, (x,y) = laurent_polynomial_ring(GF(2), [:x, :y])
            table_i = [
                (12, 4 , 1 + x + x*y      , 1 + y + x*y      , (0,  3),  (2 ,  1)),
                (14, 6 , 1 + x + y        , 1 + y + x        , (0,  7),  (1 ,  2)),
                (18, 4 , 1 + x + x*y      , 1 + y + x*y      , (0,  3),  (3 ,  0)),
                (24, 4 , 1 + x + x*y      , 1 + y + x*y      , (0,  3),  (4 ,  2)),
                (28, 6 , 1 + x + x^-1*y   , 1 + y + x*y      , (0,  7),  (2 ,  3)),
                (30, 4 , 1 + x + x^2      , 1 + y + x^2      , (0,  3),  (5 ,  1)),
                (36, 4 , 1 + x + x^-1     , 1 + y + y^-1     , (0,  9),  (2 ,  4)),
                (42, 6 , 1 + x + x*y      , 1 + y + x*y^-1   , (0,  7),  (3 ,  2)),
                (48, 4 , 1 + x + x^2      , 1 + y + x^2      , (0,  3),  (8 ,  1)),
                (54, 8 , 1 + x + x^-1     , 1 + y + x^3*y^2  , (0,  3),  (9 ,  0)),
                (56, 6 , 1 + x + y^-2     , 1 + y + x^-2     , (0,  7),  (4 ,  3)),
                (60, 8 , 1 + x + y^-2     , 1 + y + x^2      , (0, 10),  (3 ,  3)),
                (62, 10, 1 + x + x^-1*y   , 1 + y + x^-1*y^-1, (0, 31),  (1 , 13)),
                (66,  4, 1 + x + x^-2*y^-1, 1 + y + x^2*y    , (0,  3),  (11,  2)),
                (70,  6, 1 + x + x*y      , 1 + y + x*y^-1   , (0,  7),  (5 ,  1)),
                (72,  8, 1 + x + x^-1*y^3 , 1 + y + x^3*y^-1 , (0, 12),  (3 ,  3)),
                (78,  4, 1 + x + x^-2*y^-1, 1 + y + x^2*y    , (0,  3),  (13,  1)),
                (84,  6, 1 + x + x^-2     , 1 + y + x^-2*y^2 , (0, 14),  (3 , -6)),
                (90,  8, 1 + x + x^-1*y^-3, 1 + y + x^3*y^-1 , (0, 15),  (3 , -6)),
                (96,  4, 1 + x + x^-2*y   , 1 + y + x*y^-2   , (0, 12),  (4 ,  2)),
                (98,  6, 1 + x + x^-1*y^2 , 1 + y + x^-2*y^-1, (0,  7),  (7 ,  0)),
                (102, 4, 1 + x + x^-3*y   , 1 + y + x^3*y^2  , (0,  3),  (17,  2)),
                (108, 8, 1 + x + x^-1*y^-3, 1 + y + x^3*y^-1 , (0,  9),  (6 ,  0)),
                (108, 8, 1 + x + x^-1*y^3 , 1 + y + x^3*y^-1 , (0,  9),  (6 ,  0))
            ]

            table_ii = [
                (112, 6 , 1 + x + x^-1*y^2 , 1 + y + x^-2*y^-1, (0,  7),  (8 ,   2)),
                (114, 4 , 1 + x + x^-3*y   , 1 + y + x^-5     , (0,  3),  (19,   1)),
                (120, 8 , 1 + x + x^-2*y   , 1 + y + x*y^2    , (0, 10),  (6 ,   4)),
                (124, 10, 1 + x + x^-1*y^2 , 1 + y + x^-2*y^-1, (0, 31),  (2 , -12)),
                (126, 12, 1 + x + x^-1*y^-2, 1 + y + x*y^-1   , (0,  9),  (7 ,   3)),
                (132, 4 , 1 + x + y^-2     , 1 + y + x^-2     , (0, 33),  (2 ,  -7)),
                (138, 4 , 1 + x + x^-3*y   , 1 + y + x^3*y^2  , (0,  3),  (23,   2)),
                (140, 6 , 1 + x + x^-2     , 1 + y + x^-2*y^2 , (0,  7),  (10,   1)),
                (144, 12, 1 + x + x^-1*y^-3, 1 + y + x^3*y^-1 , (0, 12),  (6 ,   0)),
                (144, 12, 1 + x + x^-1*y^3 , 1 + y + x^3*y^-1 , (0, 12),  (6 ,   0)),
                (146, 18, 1 + x + y^2      , 1 + y + x^-4*y   , (0, 73),  (1 ,  16)),
                (150, 8 , 1 + x + x^-2*y   , 1 + y + x*y^2    , (0, 25),  (3 ,   7)),
                (154, 6 , 1 + x + x^-1*y^2 , 1 + y + y^-4     , (0, 77),  (1 ,  16)),
                (156, 4 , 1 + x + x^-2*y   , 1 + y + x*y^-2   , (0, 39),  (2 , -11)),
                (162, 8 , 1 + x + x^-1*y^-3, 1 + y + x^3*y^-1 , (0,  9),  (9 ,  -3)),
                (162, 8 , 1 + x + x^-1*y^3 , 1 + y + x^3*y^-1 , (0,  9),  (9 ,  -3)),
                (168, 8 , 1 + x + x^-1*y^-3, 1 + y + x^3*y^-1 , (0, 42),  (2 , -16)),
                (170, 16, 1 + x + y^-4     , 1 + y + x^4      , (0, 17),  (5 ,  -7)),
                (174, 4 , 1 + x + x^-8*y   , 1 + y + x^6*y^2  , (0,  3),  (29,   1)),
                (180, 8 , 1 + x + x^-1*y^-3, 1 + y + x^-3*y^-1, (0, 15),  (6 ,   6)),
                (180, 8 , 1 + x + x^-1*y^3 , 1 + y + x^-3*y^-1, (0, 15),  (6 ,   3)),
                (182, 6 , 1 + x + x^2*y^3  , 1 + y + x^4*y    , (0,  7),  (13,   1)),
                (186, 10, 1 + x + x^2*y^3  , 1 + y + x^2*y^-2 , (0, 31),  (3 ,   7)),
                (192, 8 , 1 + x + x^-1*y^3 , 1 + y + x^3*y^-1 , (0, 12),  (8 ,   2)),
                (196, 6 , 1 + x + x^-1*y^2 , 1 + y + x^-2*y^-1, (0, 49),  (2 ,  -10))
            ]

            table_iii = [
                (198, 8 , 1 + x + x^−4     , 1 + y + x^−3*y^2 , (0,  33), (3 ,  9)),
                (204, 4 , 1 + x + x^−3*y   , 1 + y + x^−1*y^−2, (0,  51), (2 , 14)),
                (210, 10, 1 + x + x^−3*y^2 , 1 + y + x^−3*y^−1, (0,  21), (5 , 10)),
                (216, 8 , 1 + x + x^−2*y^−5, 1 + y + x^−1*y^−3, (0,  54), (2 , 16)),
                (222, 4 , 1 + x + x^−6*y^−1, 1 + y + x^5      , (0,  3),  (37,  2)),
                (224, 6 , 1 + x + x^−3*y^2 , 1 + y + x^−3*y^−1, (0,  28), (4 , −6)),
                (228, 4 , 1 + x + x^−2*y   , 1 + y + x*y^−2   , (0,  57), (2 , 10)),
                (234, 8 , 1 + x + x^−1*y^−3, 1 + y + x^3*y^−1 , (0,  39), (3 , −9)),
                (234, 8 , 1 + x + x^−1*y^3 , 1 + y + x^3*y^−1 , (0,  39), (3 ,  6)),
                (238, 6 , 1 + x + x^−4     , 1 + y + x^−3*y^2 , (0,  7),  (17,  1)),
                (240, 8 , 1 + x + x^−2*y   , 1 + y + x*y^2    , (0,  10), (12,  3)),
                (246, 4 , 1 + x + x^3*y    , 1 + y + x^2*y^−2 , (0, 123), (1 , 22)),
                (248, 10, 1 + x + x^−2*y   , 1 + y + x^−3*y^−2, (0,  62), (2 , 25)),
                (252, 12, 1 + x + x^−3*y^−1, 1 + y + x^2*y^−2 , (0,  18), (7 ,  7)),
                (254, 14, 1 + x + x^−1*y^−3, 1 + y + y^−6     , (0, 127), (1 , 25)),
                (258, 4 , 1 + x + x^−8*y^−1, 1 + y + x^5*y    , (0,  3),  (43,  1)),
                (264, 8 , 1 + x + x*y^−5   , 1 + y + x*y^4    , (0,  66), (2 , 28)),
                (266, 6 , 1 + x + x^−1*y^−1, 1 + y + x^5      , (0,  7),  (19,  2)),
                (270, 8 , 1 + x + x^−1*y^−3, 1 + y + x^3*y^−1 , (0,  15), (9 ,  6)),
                (270, 8 , 1 + x + x^−1*y^3 , 1 + y + x^3*y^−1 , (0,  45), (3 ,−12)),
                (276, 4 , 1 + x + x^−3*y   , 1 + y + x^3*y^2  , (0,  6),  (23,  5)),
                (280, 6 , 1 + x + x*y^3    , 1 + y + x^2*y^−2 , (0,  28), (5 , 12)),
                (282, 4 , 1 + x + x^−1*y^3 , 1 + y + x^3*y^−1 , (0, 141), (1 ,  7)),
                (288, 16, 1 + x + x^−1*y^3 , 1 + y + x^3*y^−1 , (0,  12), (12,  0)),
                (288, 12, 1 + x + x^−1*y^−3, 1 + y + x^3*y^−1 , (0,  12), (12,  0)),
                (292, 18, 1 + x + y^2      , 1 + y + x^−4*y   , (0,  73), (2 , 32))
            ]

            table_iv = [
                (294, 10, 1 + x + x^−3*y   , 1 + y + x*y^−3   , (0,  21),  (7 ,   7)),
                (300, 8 , 1 + x + x^−1*y^−4, 1 + y + x^−3*y^3 , (0,  75),  (2 ,  26)),
                (306, 8 , 1 + x + x^−1*y^−3, 1 + y + x^3*y^−1 , (0,  51),  (3 ,  21)),
                (308, 6 , 1 + x + x^−1*y^−2, 1 + y + x^2*y^−1 , (0,  77),  (2 , −13)),
                (310, 10, 1 + x + x^3*y^2  , 1 + y + x^−4*y^4 , (0,  31),  (5 ,  11)),
                (312, 8 , 1 + x + x^−1*y^3 , 1 + y + x*y^3    , (0,  78),  (2 , −16)),
                (318, 4 , 1 + x + x^3*y^−4 , 1 + y + x^−1*y^−3, (0, 159),  (1 ,  17)),
                (322, 6 , 1 + x + x^−3*y^2 , 1 + y + x^−4*y^−1, (0,   7),  (23,  3)),
                (324, 8 , 1 + x + x^−1*y^−3, 1 + y + x^3*y^−1 , (0,  18),  (9 ,   6)),
                (330, 8 , 1 + x + x^−6*y^2 , 1 + y + x^2*y^5  , (0,  55),  (3 ,  23)),
                (336, 10, 1 + x + x^−4     , 1 + y + x^−1*y^−3, (0,  84),  (2 ,  37)),
                (340, 16, 1 + x + y^−4     , 1 + y + x^4      , (0,  34),  (5 ,  −7)),
                (342, 8 , 1 + x + x^−1*y^−3, 1 + y + x^3*y^−1 , (0,  57),  (3 ,  15)),
                (348, 4 , 1 + x + x^−2*y^2 , 1 + y + x^−1*y^−2, (0,  87),  (2 ,  14)),
                (350, 6 , 1 + x + x^2*y^2  , 1 + y +  x^−4*y  , (0,  35),  (5 ,  13)),
                (354, 4 , 1 + x + x^−2*y^2 , 1 + y + x^−1*y^−2, (0, 177),  (1 , −53)),
                (360, 12, 1 + x + x^−1*y^3 , 1 + y + x^3*y^−1 , (0,  30),  (6 ,   6)),
                (364, 6 , 1 + x + x^−1*y^3 , 1 + y + x^3      , (0,  14),  (13,   4)),
                (366, 4 , 1 + x + x^2*y^3  , 1 + y + x^2*y^−2 , (0, 183),  (1 ,  76)),
                (372, 10, 1 + x + x^−3*y^−2, 1 + y + x^−1*y^−3, (0,  93),  (2 , −16)),
                (378, 12, 1 + x + x^3*y^−3 , 1 + y + x^4      , (0,  21),  (9 ,   6)),
                (384, 12, 1 + x + x^−4*y^−3, 1 + y + x^3*y^−1 , (0,  48),  (4 ,  20)),
                (390, 8 , 1 + x + x^−2*y^3 , 1 + y + x^2*y^3  , (0,  15),  (13,   1)),
                (392, 6 , 1 + x + x^−3*y^2 , 1 + y + x^−3*y^−1, (0,  28),  (7 ,   7)),
                (396, 8 , 1 + x + x^−1*y^−3, 1 + y + x^3*y^−1 , (0,  66),  (3 ,  18))
            ]

            for (n, k, f, g, α1, α2) in vcat(table_i, table_ii, table_iii, table_iv)
                c = GeneralizedToricCode(f, g, α1, α2)
                stab = parity_checks(c)
                mat = matrix(GF(2), stab_to_gf2(stab))
                computed_rank = rank(mat)
                @test computed_rank == code_n(c) - code_k(c)
                @test code_n(c) == n == code_n(stab)
                @test code_k(c) == k == code_k(stab)
                @test stab_looks_good(stab, remove_redundant_rows=true) == true
            end
        end

        @testset "Appendix B: Kitaev Toric Code" begin
            R, (x, y) = laurent_polynomial_ring(GF(2), [:x, :y])
            f = 1 + x
            g = 1 + y
            α1 = (0, 6)
            α2 = (3, 3)
            n, k = 36, 2
            c = GeneralizedToricCode(f, g, α1, α2)
            stab = parity_checks(c)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == n == code_n(stab)
            @test code_k(c) == k == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
            Hx = matrix(GF(2), parity_matrix_x(c))
            Hz = matrix(GF(2), parity_matrix_z(c))
            @test rank(Hx) == 17 && rank(Hz) == 17 # B7
            @test all(sum(parity_matrix_x(c), dims=1) .== 2)
            # Each column contains exactly four ones [liang2025generalizedtoriccodestwisted](@cite)
            @test all(sum(parity_matrix_x(c), dims=2) .== 4)
            @test all(sum(parity_matrix_z(c), dims=1) .== 2)
            # Each column contains exactly four ones [liang2025generalizedtoriccodestwisted](@cite)
            @test all(sum(parity_matrix_z(c), dims=2) .== 4)
        end

        @testset "[[144, 12, 12]] Generalized Toric code properties" begin
            R, (x, y) = laurent_polynomial_ring(GF(2), [:x, :y])
            n  = 144
            k = 12
            f = 1 + x + x^-1*y^-3
            g = 1 + y + x^3*y^-1
            α1 = (0, 12)
            α2 = (6,  0)
            c = GeneralizedToricCode(f, g, α1, α2)
            stab = parity_checks(c)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == n == code_n(stab)
            @test code_k(c) == k == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
            Hx = matrix(GF(2), parity_matrix_x(c))
            Hz = matrix(GF(2), parity_matrix_z(c))
            @test all(sum(parity_matrix_x(c), dims=1) .== 3)
            @test all(sum(parity_matrix_x(c), dims=2) .== 6)
            @test all(sum(parity_matrix_z(c), dims=1) .== 3)
            @test all(sum(parity_matrix_z(c), dims=2) .== 6)
            f = 1 + x + x^-1*y^3
            g = 1 + y + x^3*y^-1
            c = GeneralizedToricCode(f, g, α1, α2)
            stab = parity_checks(c)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == n == code_n(stab)
            @test code_k(c) == k == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
            Hx = matrix(GF(2), parity_matrix_x(c))
            Hz = matrix(GF(2), parity_matrix_z(c))
            @test all(sum(parity_matrix_x(c), dims=1) .== 3)
            @test all(sum(parity_matrix_x(c), dims=2) .== 6)
            @test all(sum(parity_matrix_z(c), dims=1) .== 3)
            @test all(sum(parity_matrix_z(c), dims=2) .== 6)
        end

        @testset "Table I: Distance test for some relatively small Generalized Toric codes" begin
            R, (x, y) = laurent_polynomial_ring(GF(2), [:x, :y])
            # [[12, 4, 2]]
            f = 1 + x + x*y
            g = 1 + y + x*y
            α1 = (0, 3)
            α2 = (2, 1)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 2

            # [[14, 6, 2]]
            f = 1 + x + y
            g = 1 + y + x
            α1 = (0, 7)
            α2 = (1, 2)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 2

            # [[18, 4, 4]] 
            f = 1 + x + x*y
            g = 1 + y + x*y
            α1 = (0, 3)
            α2 = (3, 0)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        
            # [[24, 4, 4]]
            f = 1 + x + x*y
            g = 1 + y + x*y
            α1 = (0, 3)
            α2 = (4, 2)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

            # [[28, 6, 4]]
            f = 1 + x + x^-1*y
            g = 1 + y + x*y
            α1 = (0, 7)
            α2 = (2, 3)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4

            # [[30, 4, 6]]
            f = 1 + x + x^2
            g = 1 + y + x^2
            α1 = (0, 3)
            α2 = (5, 1)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

            # [[36, 4, 6]] 
            f = 1 + x + x^-1
            g = 1 + y + y^-1
            α1 = (0, 9)
            α2 = (2, 4)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

            # [[42, 6, 6]]
            f = 1 + x + x*y
            g = 1 + y + x*y^-1
            α1 = (0, 7)
            α2 = (3, 2)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

            # [[48, 4, 8]]
            f = 1 + x + x^2
            g = 1 + y + x^2
            α1 = (0, 3)
            α2 = (8, 1)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8

            # [[54, 8, 6]]
            f = 1 + x + x^-1
            g = 1 + y + x^3*y^2
            α1 = (0, 3)
            α2 = (9, 0)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

            # [[56, 6, 8]]
            f = 1 + x + y^-2
            g = 1 + y + x^-2
            α1 = (0, 7)
            α2 = (4, 3)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8

            # [[60, 8, 6]] 
            f = 1 + x + y^-2
            g = 1 + y + x^2
            α1 = (0, 10)
            α2 = (3,  3)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

            # [[62, 10, 6]]
            f = 1 + x + x^-1*y
            g = 1 + y + x^-1*y^-1
            α1 = (0, 31)
            α2 = (1, 13)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6

            # [[66, 4, 10]]
            f = 1 + x + x^-2*y^-1
            g = 1 + y + x^2*y
            α1 = (0 , 3)
            α2 = (11, 2)
            c = GeneralizedToricCode(f, g, α1, α2)
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 10
        end
    end
end
