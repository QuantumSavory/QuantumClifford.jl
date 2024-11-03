using QuantumClifford
using QuantumClifford: GeneralizedStabilizer, rowdecompose, PauliChannel, mul_left!, mul_right!
using QuantumClifford: @S_str, random_stabilizer
using QuantumOpticsBase
using LinearAlgebra
using Test
using InteractiveUtils
using Random

qo_basis = SpinBasis(1//2)
qo_tgate = sparse(identityoperator(qo_basis))
qo_tgate.data[2,2] = exp(im*pi/4)

##

@testset "expect" begin
    for s in [S"X", S"Y", S"Z", S"-X", S"-Y", S"-Z"]
        for p in [P"X", P"Y", P"Z", P"-X", P"-Y", P"-Z"]
            gs = GeneralizedStabilizer(s)
            apply!(gs, pcT)
            ρ = dm(qo_tgate*Ket(s))
            @test Operator(gs) ≈ ρ
            @test isapprox(expect(p, gs), expect(Operator(p),ρ); atol=1e-5)
        end
    end

    for _ in 1:10
        for n in 1:1
            i = rand(1:n)
            stab = random_stabilizer(n)
            genstab = GeneralizedStabilizer(stab)
            ket = Ket(stab)
            @test dm(ket) ≈ Operator(stab)
            @test dm(ket) ≈ Operator(genstab)

            pauli = random_pauli(n; nophase=false, realphase=true)
            qo_pauli = Operator(pauli)

            qo_bigtgate = n==1 ? qo_tgate : embed(qo_basis^n, i, qo_tgate)
            bigtgate = embed(n,i, pcT)
            @test qo_bigtgate ≈ Operator(bigtgate)

            for step in 1:10
                # apply!(ket, qo_bigtgate) TODO implement this API
                ket = qo_bigtgate*ket
                apply!(genstab, bigtgate)
                @test dm(ket) ≈ Operator(genstab)
                @test isapprox(expect(qo_pauli, ket), expect(pauli, genstab); atol=1e-5)
            end
        end
    end
end

@testset "Single-qubit projections of stabilizer states" begin
    checks = [(random_stabilizer(1), random_pauli(1))]
    for (s, p) in checks
        gs = GeneralizedStabilizer(s)
        apply!(gs, p)
        qo_state = Operator(gs)
        project!(gs, p)[1]
        qo_state_after_proj = Operator(gs)
        qo_pauli = Operator(p)
        qo_proj1 = (identityoperator(qo_pauli) - qo_pauli)/2
        qo_proj2 = (identityoperator(qo_pauli) + qo_pauli)/2
        result1 = qo_proj1*qo_state*qo_proj1'
        result2 = qo_proj2*qo_state*qo_proj2'
        @test result1 == zero(qo_tgate) || result2 == zero(qo_tgate) # https://github.com/QuantumSavory/QuantumClifford.jl/pull/355#discussion_r1826568296
        @test qo_state_after_proj ≈ result2 || qo_state_after_proj ≈ result1
    end
end
