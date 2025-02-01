@testitem "ECC Lacross" begin
    using LinearAlgebra
    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using QuantumClifford.ECC: Lacross, code_k, code_n, parity_checks

    @testset "Reproduce Figure 3 of [pecorari2025high](@cite)" begin
        # [[52, 4, 4]]
        n = 6
        k = 3
        coeffs = [1,1,0]
        full_rank = true
        c = parity_checks(Lacross(n,k,coeffs,full_rank))
        @test code_n(c) == 52 && code_k(c) == 4
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[65, 9, 4]]
        n = 7
        k = 3
        coeffs = [1,0,1]
        full_rank = true
        c = parity_checks(Lacross(n,k,coeffs,full_rank))
        @test code_n(c) == 65 && code_k(c) == 9
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[98, 18, 4]]
        n = 7
        k = 3
        coeffs = [1,0,1]
        full_rank = false
        c = parity_checks(Lacross(n,k,coeffs,full_rank))
        @test code_n(c) == 98 && code_k(c) == 18
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[244, 4, 8]]
        n = 12
        k = 3
        coeffs = [1,1,0]
        full_rank = true
        c = parity_checks(Lacross(n,k,coeffs,full_rank))
        @test code_n(c) == 244 && code_k(c) == 4
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[225, 9, 6]]
        n = 12
        k = 3
        coeffs = [1,1,1]
        full_rank = true
        c = parity_checks(Lacross(n,k,coeffs,full_rank))
        @test code_n(c) == 225 && code_k(c) == 9
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # cross-checks from pg.4 of https://arxiv.org/pdf/2404.13010
        # [[117, 9, d]]
        n = 9
        k = 3
        coeffs = [0,0,1]
        full_rank = true # corresponds to code with open boundary conditions
        c = parity_checks(Lacross(n,k,coeffs,full_rank))
        @test code_n(c) == 117 && code_k(c) == 9
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[162, 18, d]]
        n = 9
        k = 3
        coeffs = [0,0,1]
        full_rank = false # corresponds to code with periodic boundary conditions
        c = parity_checks(Lacross(n,k,coeffs,full_rank))
        @test code_n(c) == 162 && code_k(c) == 18
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true
    end
end
