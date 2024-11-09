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

@testset "apply!" begin
    for n in [1,2,3,4] # exponential cost in this term
        s = random_stabilizer(n)
        sm = GeneralizedStabilizer(s)
        @test dm(Ket(s)) ≈ Operator(sm)
        # test clifford gates
        for _ in 1:10
            c = random_clifford(n)
            uc = Operator(c)
            @test uc * Operator(sm) * uc' ≈ Operator(apply!(sm, c)) # sm is modified in place for the next round
        end
        # and now some (repeated) non-clifford ops
        for _ in 1:5 # exponential cost in this term
            i = rand(1:n)
            nc = embed(n, i, pcT)
            apply!(sm, nc) # in-place
            c = random_clifford(n)
            uc = Operator(c)
            @test uc * Operator(sm) * uc' ≈ Operator(apply!(sm, c)) # sm is modified in place for the next round
        end
    end
end

@testset "copy and ==" begin
    for n in 1:10
        s = random_stabilizer(n)
        sm = GeneralizedStabilizer(s)
        i = rand(1:n)
        apply!(sm, embed(n, i, pcT))
        smcopy = copy(sm)
        @test smcopy == sm
        nc = embed(n, rand(1:n), pcT)
        @test copy(nc) == nc
    end
end

@testset "Multi-qubit projections using GeneralizedStabilizer for stabilizer states" begin
    for n in 1:10
        for repetition in 1:5
            for state in [MixedDestabilizer, GeneralizedStabilizer]
                s = random_stabilizer(n)
                p = random_pauli(n)
                gs = state(s)
                apply!(gs, random_clifford(n))
                qo_state = Operator(gs)
                projectrand!(gs, p)[1]
                qo_state_after_proj = Operator(gs)
                qo_pauli = Operator(p)
                qo_proj1 = (identityoperator(qo_pauli) - qo_pauli)/2
                qo_proj2 = (identityoperator(qo_pauli) + qo_pauli)/2
                result1 = qo_proj1*qo_state*qo_proj1'
                result2 = qo_proj2*qo_state*qo_proj2'
                # Normalize to ensure consistent comparison of the projected state, independent of scaling factors
                norm_qo_state_after_proj = qo_state_after_proj/tr(qo_state_after_proj)
                norm_result1 = result1/tr(result1)
                norm_result2 = result2/tr(result2)
                @test norm_qo_state_after_proj ≈ norm_result2 || norm_qo_state_after_proj ≈ norm_result1
            end
        end
    end
end
