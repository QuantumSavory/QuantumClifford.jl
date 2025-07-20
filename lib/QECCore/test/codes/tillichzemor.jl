@testitem "Quantum Tillich-Zemor" begin

    using QuantumClifford
    using QuantumClifford: stab_looks_good
    using QuantumClifford.ECC
    using Nemo: matrix, GF
    using QECCore.LinearAlgebra
    using QECCore

    @testset "Testing Quantum Tillich-Zemor properties" begin
        for n in 4:2:20
            for m in 3:10
                for r in 3:5
                    if m ≥ r && (n - m)*r ≥ m
                        c = TillichZemor(n, m, r)
                        stab = parity_checks(c)
                        nₛ, kₛ = code_n(stab), code_k(stab)
                        H = stab_to_gf2(stab)
                        mat = matrix(GF(2), H)
                        computed_rank = rank(mat)
                        @test computed_rank == nₛ - kₛ
                        @test stab_looks_good(stab, remove_redundant_rows=true)
                        # test the derivation of [[N, K, _ ]] formules for these novel codes.
                        @test code_n(c) == nₛ && code_k(c) == kₛ
                    end
                end
            end
        end
    end
end
