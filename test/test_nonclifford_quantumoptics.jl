using QuantumClifford
using QuantumClifford: GeneralizedStabilizer, rowdecompose, PauliChannel, mul_left!, mul_right!, invsparsity
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

function _projrand(τ,p)
    qo_state = Operator(τ)
    projectrand!(τ, p)[1] # in-place
    qo_state_after_proj = Operator(τ)
    qo_pauli = Operator(p)
    qo_proj1 = (identityoperator(qo_pauli) - qo_pauli)/2
    qo_proj2 = (identityoperator(qo_pauli) + qo_pauli)/2
    result1 = qo_proj1*qo_state*qo_proj1'
    result2 = qo_proj2*qo_state*qo_proj2'
    # Normalize to ensure consistent comparison of the projected state
    norm_qo_state_after_proj = iszero(qo_state_after_proj) ? qo_state_after_proj : qo_state_after_proj/tr(qo_state_after_proj)
    norm_result1 = iszero(result1) ? result1 : result1/tr(result1)
    norm_result2 = iszero(result2) ? result2 : result2/tr(result2)
    return norm_qo_state_after_proj, norm_result1, norm_result2
end

@testset "Single-qubit projections using for stabilizer states" begin
    for s in [S"X", S"Y", S"Z", S"-X", S"-Y", S"-Z"]
        for p in [P"X", P"Y", P"Z", P"-X", P"-Y", P"-Z"]
            genstab = GeneralizedStabilizer(s)
            apply!(genstab, pcT) # in-place
            norm_qo_state_after_proj, norm_result1, norm_result2 = _projrand(genstab,p) # in-place
            !(iszero(norm_qo_state_after_proj)) && @test real(tr(norm_qo_state_after_proj)) ≈ 1
            @test projectrand!(genstab, p)[1] |> invsparsity <= genstab |> invsparsity # Λ(χ′) ≤ Λ(χ)
            @test norm_qo_state_after_proj ≈ norm_result2 || norm_qo_state_after_proj ≈ norm_result1
       end
    end

    for _ in 1:100
        for n in 1:1
            stab = random_stabilizer(n)
            genstab = GeneralizedStabilizer(stab)
            pauli = random_pauli(n)
            apply!(genstab, pcT) # in-place
            norm_qo_state_after_proj, norm_result1, norm_result2 = _projrand(genstab,pauli) # in-place
            !(iszero(norm_qo_state_after_proj)) && @test real(tr(norm_qo_state_after_proj)) ≈ 1
            @test projectrand!(genstab, pauli)[1] |> invsparsity <= genstab |> invsparsity # Λ(χ′) ≤ Λ(χ)
            @test norm_qo_state_after_proj ≈ norm_result2 || norm_qo_state_after_proj ≈ norm_result1
        end
    end
end

@testset "Multi-qubit projections using GeneralizedStabilizer for stabilizer states" begin
    for n in 1:5
        for repetition in 1:3
            for state in [Stabilizer, MixedDestabilizer, GeneralizedStabilizer]
                s = random_stabilizer(n)
                p = random_pauli(n)
                τ = state(s)
                apply!(τ, random_clifford(n)) # in-place
                norm_qo_state_after_proj, norm_result1, norm_result2 = _projrand(τ,p) # in-place
                !(iszero(norm_qo_state_after_proj)) && @test real(tr(norm_qo_state_after_proj)) ≈ 1
                @test norm_qo_state_after_proj ≈ norm_result2 || norm_qo_state_after_proj ≈ norm_result1
                isa(τ, GeneralizedStabilizer) && @test projectrand!(τ, p)[1] |> invsparsity <= τ |> invsparsity # Λ(χ′) ≤ Λ(χ)
            end
        end
    end
end

@testset "The trace Tr[χ′] is the probability of measuring an outcome" begin
    # The trace Tr[χ′] represents the probability of obtaining an outcome of 0.
    # TODO Document - Since ((real(expect(P"Z", apply!(genstab(S"-Z"), pcT))))+1)/2 is 0,
    # it triggers Eq. 16, where (I+(-1)^(i·c)*M)ρₛ(I+(-1)^(j·c)*M) evaluates to 0. Thus,
    # genstab after projectrand!(apply!(genstab(S"-Z"), pcT), P"Z")[1] has no meaning.
    for s in [S"X", S"Y", S"Z", S"-X", S"-Y", S"Z"]
        for p in [P"X", P"Y", P"Z", P"-X", P"-Y"]
            genstab = GeneralizedStabilizer(s)
            apply!(genstab, pcT) # in-place
            prob1 = (real(expect(p, genstab))+1)/2
            projectrand!(genstab, p)[1] # in-place
            dict = genstab.destabweights
            trace_χ′ = real(tr(collect(values(dict))[1])) # Tr[χ′]
            @test isapprox(prob1, trace_χ′; atol=1e-5)
        end
    end
end

@testset "Multi-qubit projections using GeneralizedStabilizer with multiple non-Clifford gates" begin
    # TODO Analyze some multi-qubit genstab states that are unsimulable due to very complex
    # destabweights, which also exhibit an inverse sparsity relation (Λ(χ′) = Λ(χ) = 4).
    count = 0
    num_trials = 5
    num_qubits = [2,3,4,5] # exclusively multi-qubit
    for n in num_qubits # exponential cost in this term
        for repetition in 1:num_trials
            stab = random_stabilizer(n)
            pauli = random_pauli(n)
            genstab = GeneralizedStabilizer(stab)
            # some (repeated) non-clifford ops
            i = rand(1:n)
            nc = embed(n, i, pcT)
            apply!(genstab, nc) # in-place
            apply!(genstab, nc) # in-place
            apply!(genstab, nc) # in-place
            norm_qo_state_after_proj, norm_result1, norm_result2 = _projrand(genstab,pauli) # in-place
            !(iszero(norm_qo_state_after_proj)) && @test real(tr(norm_qo_state_after_proj)) ≈ 1
            @test projectrand!(genstab, pauli)[1] |> invsparsity <= genstab |> invsparsity # Λ(χ′) ≤ Λ(χ)
            count += (norm_qo_state_after_proj ≈ norm_result2 || norm_qo_state_after_proj ≈ norm_result1) ? 1 : 0
        end
    end
    prob = count/(num_trials*(length(num_qubits)))
    @test prob > 0.7
end
