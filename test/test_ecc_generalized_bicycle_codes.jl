@testitem "ECC GB" tags=[:ecc] begin
    using Hecke
    using HiGHS
    using JuMP
    using QuantumClifford: stab_looks_good, stab_to_gf2
    using QuantumClifford.ECC.QECCore: code_k, code_n, distance, rate
    using QuantumClifford.ECC: generalized_bicycle_codes, code_k, code_n, DistanceMIPAlgorithm, parity_checks, GeneralizedBicycleCode, ExtendedGeneralizedBicycleCode

    # codes taken from Table 1 of [lin2024quantum](@cite)
    # Abelian 2BGA codes can be viewed as GB codes.
    @testset "GB codes" begin
        # [[70, 8, 10]]
        c = generalized_bicycle_codes([0, 15, 16, 18], [0, 1, 24, 27], 35)
        @test code_n(c) == 70 && code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_qubit=1)) == 10
        # [[54, 6, 9]]
        c = generalized_bicycle_codes([0, 1, 3, 7], [0, 1, 12, 19], 27)
        @test code_n(c) == 54 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_qubit=1)) == 9
        # [[60, 6, 10]]
        c = generalized_bicycle_codes([0 , 10, 6, 13], [0, 25, 16, 12], 30)
        @test code_n(c) == 60 && code_k(c) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_qubit=1)) == 10
        # [[72, 8, 10]]
        c = generalized_bicycle_codes([0, 9, 28, 31], [0, 1, 21, 34], 36)
        @test code_n(c) == 72 && code_k(c) == 8
        # [[72, 10, 9]]
        c = generalized_bicycle_codes([0, 9, 28, 13], [0, 1, 3, 22], 36)
        @test code_n(c) == 72 && code_k(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_qubit=1)) == 9
    end

    @testset "Codes from Appendix B of [koukoulekidis2024smallquantumcodesalgebraic](@cite))" begin
        # [[10, 2, 3]]
        R, x = polynomial_ring(GF(2), "x")
        l = 5
        a = 1 + x^4
        b = 1 + x + x^2 + x^4
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 10 == code_n(stab) && code_k(c) == 2 == code_k(stab)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test stab_looks_good(stab, remove_redundant_rows=true) == true
        for m in 2:10
            R, x = polynomial_ring(GF(2), :x)
            ext_code = ExtendedGeneralizedBicycleCode(c, m, one(R))
            stab = parity_checks(ext_code)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(ext_code) - code_k(ext_code)
            @test code_n(ext_code) == m*10 == code_n(stab) && code_k(ext_code) == 2 == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
        end

        # [[12, 2, 3]]
        R, x = polynomial_ring(GF(2), "x")
        l = 6
        a = 1 + x + x^2 + x^5
        b = 1 + x + x^3 + x^5
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 12 == code_n(stab) && code_k(c) == 2 == code_k(stab)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test stab_looks_good(stab, remove_redundant_rows=true) == true
        for m in 2:10
            R, x = polynomial_ring(GF(2), :x)
            ext_code = ExtendedGeneralizedBicycleCode(c, m, one(R))
            stab = parity_checks(ext_code)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(ext_code) - code_k(ext_code)
            @test code_n(ext_code) == m*12 == code_n(stab) && code_k(ext_code) == 2 == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
        end

        # [[14, 2, 3]]
        R, x = polynomial_ring(GF(2), "x")
        l = 7
        a = 1 + x^3
        b = 1 + x + x^3 + x^6
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 14 == code_n(stab) && code_k(c) == 2 == code_k(stab)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test stab_looks_good(stab, remove_redundant_rows=true) == true
        for m in 2:10
            R, x = polynomial_ring(GF(2), :x)
            ext_code = ExtendedGeneralizedBicycleCode(c, m, one(R))
            stab = parity_checks(ext_code)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(ext_code) - code_k(ext_code)
            @test code_n(ext_code) == m*14 == code_n(stab) && code_k(ext_code) == 2 == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
        end

        # [[16, 2, 3]]
        R, x = polynomial_ring(GF(2), "x")
        l = 8
        a = x + x^3
        b = 1 + x^5
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 16 == code_n(stab) && code_k(c) == 2 == code_k(stab)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test stab_looks_good(stab, remove_redundant_rows=true) == true
        for m in 2:10
            R, x = polynomial_ring(GF(2), :x)
            ext_code = ExtendedGeneralizedBicycleCode(c, m, one(R))
            stab = parity_checks(ext_code)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(ext_code) - code_k(ext_code)
            @test code_n(ext_code) == m*16 == code_n(stab) && code_k(ext_code) == 2 == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
        end

        # [[18, 2, 3]]
        R, x = polynomial_ring(GF(2), "x")
        l = 9
        a = 1 + x^2
        b = 1 + x^5
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 18 == code_n(stab) && code_k(c) == 2 == code_k(stab)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 3
        @test stab_looks_good(stab, remove_redundant_rows=true) == true
        for m in 2:10
            R, x = polynomial_ring(GF(2), :x)
            ext_code = ExtendedGeneralizedBicycleCode(c, m, one(R))
            stab = parity_checks(ext_code)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(ext_code) - code_k(ext_code)
            @test code_n(ext_code) == m*18 == code_n(stab) && code_k(ext_code) == 2 == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
        end

        # [[20, 2, 4]]
        R, x = polynomial_ring(GF(2), "x")
        l = 10
        a = 1 + x
        b = 1 + x^6
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 20 == code_n(stab) && code_k(c) == 2 == code_k(stab)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        @test stab_looks_good(stab, remove_redundant_rows=true) == true
        for m in 2:10
            R, x = polynomial_ring(GF(2), :x)
            ext_code = ExtendedGeneralizedBicycleCode(c, m, one(R))
            stab = parity_checks(ext_code)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(ext_code) - code_k(ext_code)
            @test code_n(ext_code) == m*20 == code_n(stab) && code_k(ext_code) == 2 == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
        end
    end

    @testset "Codes from Table I of https://arxiv.org/pdf/2502.19406" begin
        R, t = polynomial_ring(GF(2), :t)
        # [[30, 8, 4]]
        l = 15
        a = 1 + t^6 + t^13
        b = 1 + t + t^4
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 30 == code_n(stab) && code_k(c) == 8 == code_k(stab)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        @test stab_looks_good(stab, remove_redundant_rows=true) == true

        # [[62, 10, 6]]
        l = 31
        a = 1 + t + t^12
        b = 1 + t^3 + t^8
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 62 == code_n(stab) && code_k(c) == 10 == code_k(stab)
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 6
        @test stab_looks_good(stab, remove_redundant_rows=true) == true

        # [[126, 12, 10]]
        l = 63
        a = 1 + t^7 + t^8
        b = 1 + t^37 + t^43
        c = GeneralizedBicycleCode(a, b, l)
        stab = parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 126 == code_n(stab) && code_k(c) == 12 == code_k(stab)
        @test stab_looks_good(stab, remove_redundant_rows=true) == true
    end

    @testset "Generalized Bicycle Codes from Tables V-VIII of https://arxiv.org/abs/2503.03827" begin
        R, y = polynomial_ring(GF(2), :y)
        table_v = [
            (12 , 4 , 1 + y    + y^2 , 1 + y + y^2 , 6),
            (14 , 6 , 1 + y    + y^3 , 1 + y + y^3 , 7),
            (18 , 4 , 1 + y^2  + y^4 , 1 + y + y^2 , 9),
            (24 , 4 , 1 + y^2  + y^4 , 1 + y + y^2 , 12),
            (28 , 6 , 1 + y^2  + y^3 , 1 + y + y^5 , 14),
            (30 , 8 , 1 + y^2  + y^8 , 1 + y + y^4 , 15),
            (36 , 4 , 1 + y^2  + y^4 , 1 + y + y^5 , 18),
            (42 , 10, 1 + y^2  + y^10, 1 + y + y^5 , 21),
            (48 , 4 , 1 + y^5  + y^7 , 1 + y + y^5 , 24),
            (54 , 4 , 1 + y^5  + y^7 , 1 + y + y^5 , 27),
            (56 , 6 , 1 + y^3  + y^9 , 1 + y + y^5 , 28),
            (60 , 8 , 1 + y^7  + y^9 , 1 + y + y^4 , 30),
            (62 , 10, 1 + y^3  + y^8 , 1 + y + y^12, 31),
            (66 , 4 , 1 + y^2  + y^7 , 1 + y + y^11, 33),
            (70 , 6 , 1 + y^3  + y^9 , 1 + y + y^5 , 35),
            (72 , 4 , 1 + y^2  + y^7 , 1 + y + y^11, 36),
            (78 , 4 , 1 + y^5  + y^7 , 1 + y + y^8 , 39),
            (84 , 10, 1 + y^11 + y^13, 1 + y + y^5 , 42),
            (90 , 8 , 1 + y^2  + y^9 , 1 + y + y^12, 45),
            (96 , 4 , 1 + y^5  + y^7 , 1 + y + y^11, 48),
            (98 , 6 , 1 + y^4  + y^12, 1 + y + y^10, 49),
            (102, 4 , 1 + y^4  + y^8 , 1 + y + y^11, 51),
            (108, 4 , 1 + y^8  + y^10, 1 + y + y^8 , 54)
        ]

        table_vi = [
            (112, 6 , 1 + y^3  + y^15, 1 + y + y^10, 56),
            (114, 4 , 1 + y^8  + y^13, 1 + y + y^11, 57),
            (120, 8 , 1 + y^8  + y^21, 1 + y + y^12, 60),
            (124, 10, 1 + y^8  + y^11, 1 + y + y^13, 62),
            (126, 12, 1 + y^12 + y^23, 1 + y + y^8 , 63),
            (132, 4 , 1 + y^4  + y^14, 1 + y + y^14, 66),
            (138, 4 , 1 + y^8  + y^13, 1 + y + y^8 , 69),
            (140, 6 , 1 + y^10 + y^16, 1 + y + y^12, 70),
            (144, 4 , 1 + y^23 + y^28, 1 + y + y^20, 72),
            (146, 18, 1 + y^2  + y^18, 1 + y + y^9 , 73),
            (150, 8 , 1 + y^2  + y^8 , 1 + y + y^19, 75),
            (154, 6 , 1 + y^4  + y^34, 1 + y + y^19, 77),
            (156, 4 , 1 + y^11 + y^16, 1 + y + y^14, 78),
            (162, 4 , 1 + y^7  + y^11, 1 + y + y^14, 81),
            (168, 10, 1 + y^11 + y^19, 1 + y + y^17, 84),
            (170, 16, 1 + y^21 + y^25, 1 + y + y^16, 85),
            (174, 4 , 1 + y^7  + y^11, 1 + y + y^17, 87),
            (180, 8 , 1 + y^8  + y^47, 1 + y + y^34, 90),
            (182, 6 , 1 + y^9  + y^13, 1 + y + y^38, 91),
            (186, 14, 1 + y^8  + y^19, 1 + y + y^14, 93),
            (192, 4 , 1 + y^11 + y^16, 1 + y + y^14, 96),
            (196, 6 , 1 + y^12 + y^22, 1 + y + y^19, 98)
        ]

        table_vii = [
            (198, 4 , 1 + y^11 + y^16, 1 + y + y^14,  99),
            (204, 4 , 1 + y^16 + y^35, 1 + y + y^11, 102),
            (210, 14, 1 + y^11 + y^27, 1 + y + y^19, 105),
            (216, 4 , 1 + y^14 + y^22, 1 + y + y^20, 108),
            (222, 4 , 1 + y^10 + y^14, 1 + y + y^20, 111),
            (224, 6 , 1 + y^3  + y^22, 1 + y + y^31, 112),
            (228, 4 , 1 + y^7  + y^17, 1 + y + y^20, 114),
            (234, 4 , 1 + y^13 + y^29, 1 + y + y^20, 117),
            (238, 6 , 1 + y^9  + y^20, 1 + y + y^24, 119),
            (240, 8 , 1 + y^13 + y^21, 1 + y + y^19, 120),
            (246, 4 , 1 + y^13 + y^20, 1 + y + y^23, 123),
            (248, 10, 1 + y^17 + y^27, 1 + y + y^13, 124),
            (252, 12, 1 + y^25 + y^30, 1 + y + y^8 , 126),
            (254, 14, 1 + y^10 + y^37, 1 + y + y^31, 127),
            (258, 4 , 1 + y^14 + y^19, 1 + y + y^14, 129),
            (264, 4 , 1 + y^13 + y^20, 1 + y + y^17, 132),
            (266, 6 , 1 + y^12 + y^25, 1 + y + y^17, 133),
            (270, 8 , 1 + y^6  + y^23, 1 + y + y^27, 135),
            (276, 4 , 1 + y^8  + y^31, 1 + y + y^20, 138),
            (280, 6 , 1 + y^20 + y^23, 1 + y + y^17, 140),
            (282, 4 , 1 + y^10 + y^17, 1 + y + y^23, 141),
            (288, 4 , 1 + y^20 + y^25, 1 + y + y^14, 144),
            (292, 18, 1 + y^4  + y^36, 1 + y + y^9 , 146),
        ]

        table_viii = [
            (294, 10, 1 + y^19 + y^29, 1 + y + y^26, 147),
            (300, 8 , 1 + y^43 + y^52, 1 + y + y^57, 150),
            (306, 4 , 1 + y^8  + y^22, 1 + y + y^20, 153),
            (308, 6 , 1 + y^12 + y^22, 1 + y + y^26, 154),
            (310, 10, 1 + y^20 + y^43, 1 + y + y^14, 155),
            (312, 4 , 1 + y^10 + y^17, 1 + y + y^23, 156),
            (318, 4 , 1 + y^14 + y^34, 1 + y + y^38, 159),
            (322, 6 , 1 +  y^5 + y^25, 1 + y + y^24, 161),
            (324, 4 , 1 + y^11 + y^16, 1 + y + y^26, 162),
            (330, 8 , 1 + y^32 + y^38, 1 + y + y^49, 165),
            (336, 10, 1 + y^19 + y^50, 1 + y + y^5 , 168),
            (340, 16, 1 + y^4  + y^25, 1 + y + y^70, 170),
            (342, 4 , 1 + y^16 + y^23, 1 + y + y^20, 171),
            (348, 4 , 1 + y^16 + y^23, 1 + y + y^20, 174),
            (350, 6 , 1 + y^4  + y^33, 1 + y + y^24, 175),
            (354, 4 , 1 + y^19 + y^29, 1 + y + y^23, 177),
            (360, 8 , 1 + y^5  + y^25, 1 + y + y^27, 180),
            (364, 6 , 1 + y^17 + y^22, 1 + y + y^24, 182),
            (366, 4 , 1 + y^8  + y^28, 1 + y + y^26, 183),
            (372, 14, 1 + y^26 + y^34, 1 + y + y^20, 186),
            (378, 12, 1 + y^4  + y^37, 1 + y + y^25, 189),
            (384, 4 , 1 + y^16 + y^23, 1 + y + y^26, 192),
            (390, 8 , 1 + y^13 + y^37, 1 + y + y^42, 195),
            (392, 6 , 1 + y^6  + y^37, 1 + y + y^24, 196),
            (396, 4 , 1 + y^14 + y^22, 1 + y + y^32, 198)
        ]

        for (n, k, f, g, l) in vcat(table_v, table_vi, table_vii, table_viii)
            c = GeneralizedBicycleCode(f, g, l)
            stab = parity_checks(c)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == n == code_n(stab)
            @test code_k(c) == k == code_k(stab)
            @test stab_looks_good(stab, remove_redundant_rows=true) == true
        end
    end
end
