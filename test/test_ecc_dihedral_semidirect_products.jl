@testitem "ECC 2BGA Table III via semidirect products" tags=[:ecc, :ecc_bespoke_checks, :oscar_required] begin
    using Hecke: GF, group_algebra
    using QuantumClifford.ECC: code_k, code_n, two_block_group_algebra_code
    using Oscar: automorphism_group, canonical_injection, cyclic_group, describe, hom,
        normal_subgroup, order, semidirect_product, small_group_identification

    function dihedral_group_algebra(m, small_group_index)
        cyclic = cyclic_group(m)
        complement = cyclic_group(2)
        automorphisms = automorphism_group(cyclic)
        inversion = automorphisms(hom(cyclic, cyclic, [cyclic[1]], [cyclic[1]^-1]))
        action = hom(complement, automorphisms, [complement[1]], [inversion])
        group = semidirect_product(cyclic, action, complement)

        algebra = group_algebra(GF(2), group)
        r = algebra(canonical_injection(group, 1)(cyclic[1]))
        s = algebra(canonical_injection(group, 2)(complement[1]))

        @test r^m == one(r)
        @test s^2 == one(s)
        @test (r * s)^2 == one(r)
        @test normal_subgroup(group) == cyclic
        @test order(group) == 2 * m
        @test describe(group) == "D$(2 * m)"
        @test small_group_identification(group) == (2 * m, small_group_index)

        return r, s
    end

    function check_table_entry(a_elements, b_elements, expected)
        @test length(a_elements) == 2
        @test allunique(a_elements)
        @test length(b_elements) == 6
        @test allunique(b_elements)

        code = two_block_group_algebra_code(sum(a_elements), sum(b_elements))
        @test (code_n(code), code_k(code)) == expected
    end

    @testset "Construct Table III of arXiv:2306.16400" begin
        r, s = dihedral_group_algebra(6, 4)
        # [[24, 8, 3]]
        check_table_entry(
            [one(r), r^4],
            [one(r), s * r^4, r^3, r^4, s * r^2, r],
            (24, 8),
        )
        # [[24, 12, 2]]
        check_table_entry(
            [one(r), r^3],
            [one(r), s * r, r^3, r^4, s * r^4, r],
            (24, 12),
        )

        r, s = dihedral_group_algebra(8, 7)
        # [[32, 8, 4]]
        check_table_entry(
            [one(r), r^2],
            [one(r), s * r^5, s * r^4, r^2, s * r^7, s * r^6],
            (32, 8),
        )
        # [[32, 16, 2]]
        check_table_entry(
            [one(r), r^4],
            [one(r), s * r^3, s * r^6, r^4, s * r^7, s * r^2],
            (32, 16),
        )

        r, s = dihedral_group_algebra(9, 1)
        # [[36, 12, 3]]
        check_table_entry(
            [one(r), r^3],
            [one(r), s, r, r^3, s * r^3, r^4],
            (36, 12),
        )

        r, s = dihedral_group_algebra(10, 4)
        # [[40, 8, 5]]
        check_table_entry(
            [one(r), r^2],
            [one(r), s * r^4, r^5, r^2, s * r^6, r],
            (40, 8),
        )
        # [[40, 20, 2]]
        check_table_entry(
            [one(r), r^5],
            [one(r), s * r^2, r^5, r^6, s * r^7, r],
            (40, 20),
        )

        r, s = dihedral_group_algebra(12, 6)
        # [[48, 8, 6]]
        check_table_entry(
            [one(r), r^10],
            [one(r), s * r^8, r^9, r^4, s * r^2, r^5],
            (48, 8),
        )
        # [[48, 12, 4]]
        check_table_entry(
            [one(r), r^3],
            [one(r), s * r^7, r^3, r^4, s * r^10, r^7],
            (48, 12),
        )
        # [[48, 16, 3]]
        check_table_entry(
            [one(r), r^8],
            [one(r), s * r^8, r^9, r^8, s * r^4, r^5],
            (48, 16),
        )
        # [[48, 24, 2]]
        check_table_entry(
            [one(r), r^6],
            [one(r), s * r^11, r^6, s * r^5, r, r^7],
            (48, 24),
        )

        r, s = dihedral_group_algebra(14, 3)
        # [[56, 8, 7]]
        check_table_entry(
            [one(r), r^4],
            [one(r), s * r^11, r^7, s * r^5, r^12, r^9],
            (56, 8),
        )
        # [[56, 28, 2]]
        check_table_entry(
            [one(r), r^7],
            [one(r), s * r^2, r^7, r^8, s * r^9, r],
            (56, 28),
        )

        r, s = dihedral_group_algebra(15, 3)
        # [[60, 12, 5]]
        check_table_entry(
            [one(r), r^12],
            [one(r), s * r^14, r^5, r^12, s * r^11, r^14],
            (60, 12),
        )
        # [[60, 20, 3]]
        check_table_entry(
            [one(r), r^5],
            [one(r), s * r^13, r^5, r^12, s * r^3, r^2],
            (60, 20),
        )

        r, s = dihedral_group_algebra(16, 18)
        # [[64, 8, 8]]
        check_table_entry(
            [one(r), r^6],
            [one(r), s * r^12, s * r^9, r^6, s, s * r],
            (64, 8),
        )
        # [[64, 16, 8]]
        check_table_entry(
            [one(r), r^4],
            [one(r), s * r^10, s * r^3, r^4, s * r^14, s * r^7],
            (64, 16),
        )
        # [[64, 32, 2]]
        check_table_entry(
            [one(r), r^8],
            [one(r), s * r^11, s * r^12, r^8, s * r^3, s * r^4],
            (64, 32),
        )
    end
end
