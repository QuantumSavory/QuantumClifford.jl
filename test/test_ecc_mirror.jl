@testitem "ECC Mirror Code" tags=[:ecc, :ecc_bespoke_checks, :oscar_required] begin
    using Oscar
    using QECCore
    using QuantumClifford: Stabilizer, canonicalize!, phases, stab_looks_good
    using QuantumClifford.ECC

    function contains_negative_identity(checks)
        canonical, _, rank = canonicalize!(copy(checks), ranks=true)
        any(==(0x02), @view phases(canonical)[rank+1:end])
    end

    @testset "Published Abelian mirror codes" begin
        cases = [
            # [[16, 4, 4]]
            (16, 4, [16], [(0,), (4,)], [(1,), (3,), (5,), (11,)]),
            # [[18, 4, 4]]
            (18, 4, [2, 9], [(0, 0), (0, 1), (0, 2)], [(1, 0), (1, 1), (1, 5)]),
            # [[36, 6, 6]]
            (
                36,
                6,
                [2, 2, 3, 3],
                [(0, 0, 0, 0), (0, 1, 0, 1), (1, 0, 0, 2)],
                [(0, 0, 0, 0), (0, 1, 1, 0), (1, 1, 2, 0)],
            ),
        ]

        for (n, k, orders, A, B) in cases
            G = abelian_group(orders)
            code = Mirror(G, A, B)
            H = parity_matrix(code)

            @test eltype(H) === Bool
            @test size(H) == (n, 2n)
            @test code_n(code) == n
            @test code_k(code) == k
            @test maximum(count, eachrow(H)) <= length(code.A) + length(code.B)
            @test stab_looks_good(parity_checks(code); remove_redundant_rows=true)
            @test !contains_negative_identity(parity_checks(code))

            asymmetric = Mirror(G, A, B; symmetric=false)
            @test parity_matrix(asymmetric) == H
        end
    end

    @testset "Non-Abelian checks and phases" begin
        G = dihedral_group(8)
        r = only(filter(g -> order(g) == 4, gens(G)))
        s = first(filter(g -> order(g) == 2, gens(G)))

        # This valid code has negative-identity dependencies if every row is
        # assigned a positive sign.
        code = Mirror(G, [one(G), s], [one(G), r^3, s, r*s])
        all_positive = Stabilizer(parity_matrix(code))
        corrected = parity_checks(code)
        @test contains_negative_identity(all_positive)
        @test any(==(0x02), phases(corrected))
        @test !contains_negative_identity(corrected)
        @test stab_looks_good(corrected; remove_redundant_rows=true)

        symmetric = Mirror(G, [r], [r])
        asymmetric = Mirror(G, [r], [r]; symmetric=false)
        @test parity_matrix(symmetric) != parity_matrix(asymmetric)
        @test stab_looks_good(parity_checks(asymmetric); remove_redundant_rows=true)

        @test_throws ArgumentError Mirror(G, [one(G)], [s])
        @test_throws ArgumentError Mirror(G, [one(G)], [s], true)
    end

    @testset "Input forms and validation" begin
        G = abelian_group([4])
        elements = collect(G)
        from_elements = Mirror(G, elements[1:2], elements[2:3])
        @test all(parent(element) === G for element in from_elements.A)

        deduplicated = Mirror(G, [(0,), (0,)], [(1,), (1,)])
        @test length(deduplicated.A) == length(deduplicated.B) == 1

        empty_supports = Mirror(G, [], [])
        @test code_k(empty_supports) == code_n(empty_supports) == 4

        other_G = abelian_group([4])
        @test_throws ArgumentError Mirror(G, [first(collect(other_G))], [])
        @test_throws ArgumentError Mirror(abelian_group([2, 3]), [(0,)], [])
        infinite_G = free_group(1)
        @test_throws ArgumentError Mirror(infinite_G, [one(infinite_G)], [])

        multiplicative_G = cyclic_group(4)
        g = first(gens(multiplicative_G))
        multiplicative_code = Mirror(multiplicative_G, [one(multiplicative_G), g], [g])
        @test size(parity_matrix(multiplicative_code)) == (4, 8)
    end
end
