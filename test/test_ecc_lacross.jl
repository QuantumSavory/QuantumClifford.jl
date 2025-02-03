@testitem "ECC Lacross" begin
    using LinearAlgebra
    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using QuantumClifford.ECC: Lacross, code_k, code_n, parity_checks

    @testset "Reproduce Figure 3 of [pecorari2025high](@cite)" begin
        @testset "Reproduce Figure 3(a) of [pecorari2025high](@cite)" begin
            # [[52, 4, 4]]
            n = 6
            k = 2
            coeffs = [1,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 52 && code_k(c) == 4
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[100, 4, 5]]
            n = 8
            k = 2
            coeffs = [0,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 100 && code_k(c) == 4
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[130, 4, 6]]
            n = 9
            k = 2
            coeffs = [1,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 130 && code_k(c) == 4
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[164, 4, 6]]
            n = 10
            k = 2
            coeffs = [0,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 164 && code_k(c) == 4
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[244, 4, 8]]
            n = 12
            k = 2
            coeffs = [1,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 244 && code_k(c) == 4
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true
        end

        @testset "Reproduce Figure 3(b) of [pecorari2025high](@cite)" begin
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

            # [[117, 9, 4]]
            n = 9
            k = 3
            coeffs = [0,0,1]
            full_rank =  true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 117 && code_k(c) == 9
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

            # [[317, 9, 8]]
            n = 14
            k = 3
            coeffs = [1,0,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
             @test code_n(c) == 317 && code_k(c) == 9
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[369, 9, 8]]
            n = 15
            k = 3
            coeffs = [0,0,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
             @test code_n(c) == 369 && code_k(c) == 9
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[52, 4, 4]]
            n = 6
            k = 3
            coeffs = [1,1,0]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 52 && code_k(c) == 4
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true
        end

        @testset "Reproduce Figure 3(c) of [pecorari2025high](@cite)" begin
            # [[136, 16, 5]]
            n = 10
            k = 4
            coeffs = [1,1,1,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 136 && code_k(c) == 16
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[208, 16, 6]]
            n = 12
            k = 4
            coeffs = [0,0,0,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 208 && code_k(c) == 16
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[296, 16, 7]]
            n = 14
            k = 4
            coeffs = [1,1,0,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 296 && code_k(c) == 16
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

            # [[400, 16, 8]]
            n = 16
            k = 4
            coeffs = [0,0,0,1]
            full_rank = true
            c = parity_checks(Lacross(n,k,coeffs,full_rank))
            @test code_n(c) == 400 && code_k(c) == 16
            @test stab_looks_good(copy(c), remove_redundant_rows=true) == true
        end

        @testset "cross-checks from pg.4 of https://arxiv.org/pdf/2404.13010" begin
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
end
