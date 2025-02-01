@testitem "ECC Lacross" begin
    using Hecke
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC: Lacross code_k, code_n

    @testset "Reproduce Figure 3 of [pecorari2025high](@cite)" begin
        # [[52, 4, 4]]
        n = 6
        pattern = [1,1,1,0]
        full_rank = true
        c = parity_checks(Lacross(n,pattern,full_rank))
        @test code_n(c) == 52 && code_k(c) == 4
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[65, 9, 4]]
        n = 7
        pattern = [1,1,0,1]
        full_rank = true
        c = parity_checks(Lacross(n,pattern,full_rank))
        @test code_n(c) == 65 && code_k(c) == 9
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[98, 18, 4]]
        n = 7
        pattern = [1,1,0,1]
        full_rank = false
        c = parity_checks(Lacross(n,pattern,full_rank))
        @test code_n(c) == 98 && code_k(c) == 18
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[244,4,8]]
        n = 12
        pattern = [1,1,1,0]
        full_rank = true
        c = parity_checks(Lacross(n,pattern,full_rank))
        @test code_n(c) == 244 && code_k(c) == 4
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true

        # [[225, 9, 6]]
        n = 12
        pattern = [1,1,1,1]
        full_rank = true
        c = parity_checks(Lacross(n,pattern,full_rank))
        @test code_n(c) == 225 && code_k(c) == 9
        @test stab_looks_good(copy(c), remove_redundant_rows=true) == true
end
