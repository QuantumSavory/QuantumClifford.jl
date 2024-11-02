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

@testset "Conjugate destabs" begin
    test_cases = [
        (1, [tHadamard, tPhase, tId1]),
        (2, [tCNOT, tCPHASE, tSWAP]),
        (3, [enumerate_cliffords(3, clifford_cardinality(3)), CliffordOperator(sHadamard(3), 3), CliffordOperator(sCNOT(1, 2), 3)]),
        (4, [enumerate_cliffords(4, clifford_cardinality(4)), CliffordOperator(sHadamard(4), 4), CliffordOperator(sCNOT(2, 1), 4)]),
        (5, [enumerate_cliffords(5, clifford_cardinality(5)), CliffordOperator(sHadamard(5), 5), CliffordOperator(sCNOT(2, 3), 5)])
    ]

    for (num_qubits, gates) in test_cases
        @testset "Conjugate destabs test using $num_qubits-qubit Clifford gate" begin
            for _ in 1:10
                s = random_stabilizer(num_qubits)
                sm = GeneralizedStabilizer(s)
                @test dm(Ket(s)) ≈ Operator(sm)
                for C in gates
                    @test Operator(C) * Operator(sm) * Operator(C)' ≈ Operator(apply!(sm, C))
                end
            end
        end
    end

    @testset "Conjugate destabs test using non-Clifford gate" begin
        for n in 5:10
            i = rand(1:(n-1))
            eg = embed(n, i, pcT)
            sm = GeneralizedStabilizer(random_stabilizer(n))
            @test dm(Ket(sm.stab)) ≈ Operator(sm)
            @test Operator(eg) * Operator(sm) * Operator(eg)' ≈ Operator(apply!(sm, eg))
        end
    end
end
