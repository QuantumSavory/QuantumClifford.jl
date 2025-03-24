@testitem "ECC [[8rp, (8r − 2)p − 2m, 4]] Delfosse-Reichardt code" begin

    using LinearAlgebra
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: DelfosseReichardt, _generalize_delfosse_reichardt_code
    using Nemo: matrix, GF, echelon_form

    @testset "Testing [[8rp, (8r − 2)p − 2m, 4]] DelfosseRepCode properties" begin
        # TODO use MIP solver to test minimum distance
        @testset "test [16p, 14p − 8, 4]] code family that uses RM(2,4)" begin
            r = 2
            m = 4
            for i in 2:100
                p = i
                n = 8*r*p
                k = (8*r-2)*p-2*m
                stab = parity_checks(DelfosseReichardt(p,r,m))
                H = stab_to_gf2(stab)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k
            end
        end

        # TODO use MIP solver to test minimum distance
        @testset "test [[8p, 6(p−1), 4]] code family that uses RM(1,3)" begin
            r = 1
            m = 3
            for i in 2:100
                p = i
                n = 8*r*p
                k = (8*r-2)*p-2*m
                stab = parity_checks(DelfosseReichardt(p,r,m))
                H = stab_to_gf2(stab)
                mat = matrix(GF(2), H)
                computed_rank = rank(mat)
                @test computed_rank == n - k
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
