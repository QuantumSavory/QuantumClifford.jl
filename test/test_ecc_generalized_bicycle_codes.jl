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
end
