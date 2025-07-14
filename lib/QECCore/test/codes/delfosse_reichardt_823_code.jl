@testitem "[[8p, 4p − 2, 3]] Delfosse-Reichardt Generalized [[8,2,3]] codes" begin

    using JuMP
    using HiGHS
    using QuantumClifford
    using QuantumClifford: stab_looks_good, Stabilizer
    using QuantumClifford.ECC
    using Nemo: matrix, GF
    using QECCore.LinearAlgebra
    using QECCore

    function _consistency_check(p)
        H = parity_checks(SmallestColorCode())
        rows, cols = size(H)
        tab = zero(Stabilizer, rows - 2, cols)
        H_rep₁ = parity_checks(SmallestColorCode())[1:4, :]
        H_rep₂ = parity_checks(SmallestColorCode())[5:6, :]
        rows = [hcat(fill(tab, i - 1)..., H_rep₁, fill(tab, p - i)...) for i in 1:p]
        D = vcat(rows...)
        E = hcat(fill(H_rep₂, p)...)
        extended_H = vcat(D, E)
        return extended_H
    end

    @testset "Testing [[8p, 4p − 2, 3]] Delfosse-Reichardt Generalized [[8,2,3]] code properties" begin
        for p in 1:10
            c = DelfosseReichardt823(p)
            n, k = code_n(c), code_k(c)
            stab = parity_checks(c)
            nₛ, kₛ = code_n(stab), code_k(stab)
            H = stab_to_gf2(stab)
            mat = matrix(GF(2), H)
            computed_rank = rank(mat)
            @test stab == _consistency_check(p)
            @test computed_rank == n - k && computed_rank == nₛ - kₛ && n == nₛ && k == kₛ
            @test stab_looks_good(stab, remove_redundant_rows=true)
        end
    end
end
