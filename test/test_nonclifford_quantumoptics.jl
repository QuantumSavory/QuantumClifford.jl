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
