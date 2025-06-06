@testitem "Quantum Reed-Muller" begin
    using Test
    using Nemo: echelon_form, matrix, GF
    using LinearAlgebra
    using QuantumClifford
    using QuantumClifford: canonicalize!, Stabilizer, stab_to_gf2
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC, QuantumReedMuller, Steane7, CSS
    using QuantumClifford.ECC.QECCore: code_k, code_n, distance, rate

    function designed_distance(mat)
        dist = 3
        for row in eachrow(mat)
            count = sum(row)
            if count < dist
                return false
            end
        end
        return true
    end

    @testset "Test QuantumReedMuller(r,m) properties" begin
        for m in 3:10
            stab = parity_checks(QuantumReedMuller(m))
            H = stab_to_gf2(stab)
            @test designed_distance(H) == true
            # QuantumReedMuller(3) is the Steane7 code.
            @test canonicalize!(parity_checks(Steane7())) == parity_checks(QuantumReedMuller(3))
            @test code_n(QuantumReedMuller(m)) == 2^m - 1
            @test code_k(QuantumReedMuller(m)) == 1
            @test distance(QuantumReedMuller(m)) == 3
            @test H == stab_to_gf2(parity_checks(CSS(parity_checks_x(QuantumReedMuller(m)), parity_checks_z(QuantumReedMuller(m)))))
            # [[15,1,3]] qrm code from table 1 of https://arxiv.org/pdf/1705.0010
            qrm₁₅₁₃ = S"ZIZIZIZIZIZIZIZ
                        IZZIIZZIIZZIIZZ
                        IIIZZZZIIIIZZZZ
                        IIIIIIIZZZZZZZZ
                        IIZIIIZIIIZIIIZ
                        IIIIZIZIIIIIZIZ
                        IIIIIZZIIIIIIZZ
                        IIIIIIIIIZZIIZZ
                        IIIIIIIIIIIZZZZ
                        IIIIIIIIZIZIZIZ
                        XIXIXIXIXIXIXIX
                        IXXIIXXIIXXIIXX
                        IIIXXXXIIIIXXXX
                        IIIIIIIXXXXXXXX"
            @test canonicalize!(parity_checks(qrm₁₅₁₃)) == canonicalize!(parity_checks(QuantumReedMuller(4)))
        end
    end
end
