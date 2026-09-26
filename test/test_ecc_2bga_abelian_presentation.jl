@testitem "ECC 2BGA Table II via abelian group presentations" tags=[:ecc, :ecc_bespoke_checks, :oscar_required] begin
    using HiGHS
    using JuMP
    using Hecke: GF, gens, group_algebra, quo
    using QuantumClifford.ECC: DistanceMIPAlgorithm, code_k, code_n, distance, twobga_from_fp_group
    using Oscar: describe, free_group, order, small_group_identification

    function abelian_group_presentation(m, small_group_index)
        F = free_group(["x", "s"])
        x, s = gens(F)
        G, = quo(F, [x^m, s^2, x * s * x^-1 * s^-1])

        @test order(G) == 2 * m
        @test describe(G) == "C$m x C2"
        @test small_group_identification(G) == (2 * m, small_group_index)

        x, s = gens(G)
        return G, group_algebra(GF(2), G), x, s
    end

    function check_table_entry(algebra, a_elements, b_elements, expected_parameters)
        @test length(a_elements) == 2
        @test allunique(a_elements)
        @test length(b_elements) == 6
        @test allunique(b_elements)

        code = twobga_from_fp_group(a_elements, b_elements, algebra)
        @test (code_n(code), code_k(code)) == expected_parameters
        return code
    end

    @testset "Construct Table II of arXiv:2306.16400" begin
        G, algebra, x, s = abelian_group_presentation(4, 2)
        code = check_table_entry(
            algebra,
            [one(G), x],
            [one(G), x, s, x^2, s * x, s * x^3],
            (16, 2),
        )
        # The direct Table II suite checks the distances of all 17 entries.
        @test distance(code, DistanceMIPAlgorithm(solver=HiGHS)) == 4
        check_table_entry(
            algebra,
            [one(G), x],
            [one(G), x, s, x^2, s * x, x^3],
            (16, 4),
        )
        check_table_entry(
            algebra,
            [one(G), s],
            [one(G), x, s, x^2, s * x, s * x^2],
            (16, 8),
        )

        G, algebra, x, s = abelian_group_presentation(6, 5)
        check_table_entry(
            algebra,
            [one(G), x],
            [one(G), x^3, s, x^4, x^2, s * x],
            (24, 4),
        )
        check_table_entry(
            algebra,
            [one(G), x^3],
            [one(G), x^3, s, x^4, s * x^3, x],
            (24, 12),
        )

        G, algebra, x, s = abelian_group_presentation(8, 5)
        check_table_entry(
            algebra,
            [one(G), x^6],
            [one(G), s * x^7, s * x^4, x^6, s * x^5, s * x^2],
            (32, 8),
        )
        check_table_entry(
            algebra,
            [one(G), s * x^4],
            [one(G), s * x^7, s * x^4, x^6, x^3, s * x^2],
            (32, 16),
        )

        G, algebra, x, s = abelian_group_presentation(10, 5)
        check_table_entry(
            algebra,
            [one(G), x],
            [one(G), x^5, x^6, s * x^6, x^7, s * x^3],
            (40, 4),
        )
        check_table_entry(
            algebra,
            [one(G), x^6],
            [one(G), x^5, s, x^6, x, s * x^2],
            (40, 8),
        )
        check_table_entry(
            algebra,
            [one(G), x^5],
            [one(G), x^5, s, x^6, s * x^5, x],
            (40, 20),
        )

        G, algebra, x, s = abelian_group_presentation(12, 9)
        check_table_entry(
            algebra,
            [one(G), s * x^10],
            [one(G), x^3, s * x^6, x^4, x^7, x^8],
            (48, 8),
        )
        check_table_entry(
            algebra,
            [one(G), x^3],
            [one(G), x^3, s * x^6, x^4, s * x^9, x^7],
            (48, 12),
        )
        check_table_entry(
            algebra,
            [one(G), x^4],
            [one(G), x^3, s * x^6, x^4, x^7, s * x^10],
            (48, 16),
        )
        check_table_entry(
            algebra,
            [one(G), s * x^6],
            [one(G), x^3, s * x^6, x^4, s * x^9, s * x^10],
            (48, 24),
        )

        G, algebra, x, s = abelian_group_presentation(14, 4)
        check_table_entry(
            algebra,
            [one(G), x],
            [one(G), x^7, s * x^8, x^2, x^3, s * x^11],
            (56, 4),
        )
        check_table_entry(
            algebra,
            [one(G), x^8],
            [one(G), x^7, s, x^8, x^9, s * x^4],
            (56, 8),
        )
        check_table_entry(
            algebra,
            [one(G), x^7],
            [one(G), x^7, s, x^8, s * x^7, x],
            (56, 28),
        )
    end
end
