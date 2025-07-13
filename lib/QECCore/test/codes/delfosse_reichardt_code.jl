@testitem "[[8rp, (8r − 2)p − 2m, 4]] Delfosse-Reichardt code" begin

    using JuMP
    using HiGHS
    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using Nemo: matrix, GF, echelon_form
    using QECCore
    using QECCore.LinearAlgebra
    using QECCore: _generalize_delfosse_reichardt_code, search_self_orthogonal_rm_codes

    @testset "Testing Delfosse-Reichardt code properties" begin
        # from https://arxiv.org/pdf/2008.05051"
        code_families = [
            (2, 4),  # [[16p, 14p-8, 4]] code family
            (1, 3),  # [[8p, 6(p-1), 4]] code family
            (2, 5),  # not from paper, but arbitrary to demonstrate generalization
        ]
        for (r, m) in code_families
            @testset "DelfosseReichardt with RM(r=$r, m=$m) seed" begin
                for p in 2:5
                    c = DelfosseReichardt(p, r, m)
                    stab = parity_checks(c)
                    nₛ, kₛ = code_n(stab), code_k(stab)
                    H = stab_to_gf2(stab)
                    mat = matrix(GF(2), H)
                    computed_rank = rank(mat)
                    @test computed_rank == nₛ - kₛ
                    @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
                    @test stab_looks_good(stab, remove_redundant_rows=true)
                end
             end
        end

        max_m = 10
        parameters = search_self_orthogonal_rm_codes(max_m)
        for (r, m) in parameters
            @testset "Delfosse-Reichardt code with RM(r=$r, m=$m) seed" begin
                for p in 2:5
                    code = DelfosseReichardt(p, r, m)
                    stab = parity_checks(code)
                    nₛ, kₛ = code_n(stab), code_k(stab)
                    H = stab_to_gf2(stab)
                    mat = matrix(GF(2), H)
                    computed_rank = rank(mat)
                    @test computed_rank == nₛ - kₛ
                    @test stab_looks_good(stab, remove_redundant_rows=true)
                end
            end
        end

        # crosscheck parity check matrices from VIII. 2 Generalizing the [8, 4, 4]
        # and [16, 11, 4] classical codes, pg. 11 of https://arxiv.org/pdf/2008.05051
        r = 1
        m = 3
        blocks = 2
        mat_paper = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
                     0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
                     0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
                     0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1;
                     0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
       @test echelon_form(matrix(GF(2), Matrix{Int}(_generalize_delfosse_reichardt_code(blocks,r,m)))) == echelon_form(matrix(GF(2), mat_paper))
    end
end
