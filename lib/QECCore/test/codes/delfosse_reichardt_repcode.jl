@testitem "[[4p, 2(p − 2), 4]] Delfosse-Reichardt RepCode" begin

    using JuMP
    using HiGHS
    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using Nemo: matrix, GF
    using QECCore.LinearAlgebra
    using QECCore
    using QECCore: _extend_414_repetition_code

    @testset "Testing [[4p, 2(p − 2), 4]] Delfosse-Reichardt RepCode properties" begin
        # For p = 2, the code is trivial with k = 0 so we skip it in testing since the distance of a trivial code is infinite.
        for i in 2:5
            p = 2*i
            n = 4*p
            k = 2*(p - 2)
            c = DelfosseReichardtRepCode(p)
            stab = parity_checks(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            H = stab_to_gf2(stab)
            mat = matrix(GF(2), H)
            computed_rank = rank(mat)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 4
            @test stab_looks_good(stab, remove_redundant_rows=true)
        end

        # crosscheck parity check matrices from VIII. 1 Generalizing the [4, 1, 4]
        # classical repetition code, pg. 10 of https://arxiv.org/pdf/2008.05051
        blocks = 1
        @test _extend_414_repetition_code(blocks) == [1  1  1  1;
                                                      0  0  1  1;
                                                      0  1  0  1];
        blocks = 4
        # [16, 10, 4] classical linear code
        @test _extend_414_repetition_code(blocks) == [1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0;
                                                      0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0;
                                                      0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0;
                                                      0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1;
                                                      0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1;
                                                      0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1]
        blocks = 6
        # [24, 16, 4] classical linear code
        @test _extend_414_repetition_code(blocks) == [1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                                                      0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                                                      0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0;
                                                      0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0;
                                                      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0;
                                                      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1;
                                                      0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1;
                                                      0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1];
    end
end
