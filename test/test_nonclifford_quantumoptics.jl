@testitem "Non-Clifford Quantum Optics" begin
    using QuantumClifford
    using QuantumClifford: GeneralizedStabilizer, rowdecompose, PauliChannel, mul_left!, mul_right!, invsparsity, _projectrand_notnorm, mixed_destab_looks_good, tr
    using QuantumClifford: @S_str, random_stabilizer
    using QuantumOpticsBase
    using LinearAlgebra
    using LinearAlgebra: tr
    using Test
    using InteractiveUtils
    using Random

    qo_basis = SpinBasis(1//2)
    qo_tgate = sparse(identityoperator(qo_basis))
    qo_tgate.data[2,2] = exp(im*pi/4)

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
                @test invsparsity(projectrand!(genstab, p)[1]) <= invsparsity(genstab) # Λ(χ′) ≤ Λ(χ)
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
                    isa(τ, GeneralizedStabilizer) && @test invsparsity(projectrand!(τ, p)[1]) <= invsparsity(τ) # Λ(χ′) ≤ Λ(χ)
                end
            end
        end
    end

    @testset "The trace Tr[χ′] is the probability of measuring 0 wrt _projectrand_notnorm." begin
        for s in [S"X", S"Y", S"Z", S"-X", S"-Y", S"-Z"]
            for p in [P"X", P"Y", P"Z", P"-X", P"-Y", P"-Z"]
                genstab = GeneralizedStabilizer(s)
                apply!(genstab, pcT)
                # Calculate theoretical probability
                prob1 = (real(expect(p, genstab)) + 1) / 2
                # Perform unnormalized projection and get the updated state
                genstab_proj, _ = _projectrand_notnorm(genstab, p, 0)
                if !isempty(genstab_proj.destabweights)
                    # Calculate trace of χ′
                    trace_χ′ = tr(genstab_proj)
                    # Verify Tr[χ′] equals measurement probability using _projectrand_notnorm
                    @test isapprox(prob1, real(trace_χ′); atol=1e-5)
                else
                    # Verify zero probability case
                    @test isapprox(prob1, 0.0; atol=1e-5)
                end
            end
        end
        # Test non-trivial states
        for n in 2:5
            for trial in 1:10
                s = random_stabilizer(n)
                genstab = GeneralizedStabilizer(s)
                p = random_pauli(n)
                # Apply non-Clifford operations
                for _ in 1:n
                    apply!(genstab, embed(n, rand(1:n), pcT))
                end
                # Calculate theoretical probability
                prob1 = (real(expect(p, genstab)) + 1) / 2
                # Perform unnormalized projection and get the updated state
                genstab_proj, _ = _projectrand_notnorm(genstab, p, 0)
                if !isempty(genstab_proj.destabweights)
                    # Calculate trace of χ′
                    trace_χ′ = tr(genstab_proj)
                    # Verify Tr[χ′] equals measurement probability using _projectrand_notnorm
                    @test isapprox(prob1, real(trace_χ′); atol=1e-5)
                else
                    # Verify zero probability case
                    @test isapprox(prob1, 0.0; atol=1e-5)
                end
            end
        end
    end

    @testset "Multi-qubit projections using GeneralizedStabilizer with multiple non-Clifford gates" begin
        num_trials = 10
        num_qubits = [2,3,4,5] # exclusively multi-qubit
        for n in num_qubits  # Exponential cost in this term
            for repetition in 1:num_trials
                stab = random_stabilizer(n)
                pauli = random_pauli(n)
                genstab = GeneralizedStabilizer(stab)
                # Apply some (repeated) non-Clifford operations
                i = rand(1:n)
                nc = embed(n, i, pcT)
                apply!(genstab, nc) # in-place
                apply!(genstab, nc) # in-place
                apply!(genstab, nc) # in-place
                norm_qo_state_after_proj, norm_result1, norm_result2 = _projrand(genstab, pauli)
                !(iszero(norm_qo_state_after_proj)) && @test real(tr(norm_qo_state_after_proj)) ≈ 1
                @test norm_qo_state_after_proj ≈ norm_result2 || norm_qo_state_after_proj ≈ norm_result1
            end
        end
        for j in 1:10
            num_qubits = [2,3,4,5] # exclusively multi-qubit
            for n in num_qubits # exponential cost in this term
                genstab = GeneralizedStabilizer(random_stabilizer(n))
                p = random_pauli(n)
                for i in 1:n
                    apply!(genstab, embed(n, rand(1:n), pcT))
                end
                projectrand!(genstab, p)
                # Check the trace after normalization
                trace = tr(genstab)
                !iszero(trace) && @assert trace ≈ 1
            end
        end
    end

    @testset "projectrand! gives random results" begin
        s1 = projectrand!(GeneralizedStabilizer(S"X"), P"Z")[1].stab
        @test any(projectrand!(GeneralizedStabilizer(S"X"), P"Z")[1].stab != s1 for _ in 1:10)
    end

    @testset "Tensor products of generalized stabilizers" begin
        num_trials = 3
        num_qubits = [2,3] # exclusively multi-qubit
        for n in num_qubits
            for repetition in 1:num_trials
                stab1 = random_stabilizer(n)
                genstab1 = GeneralizedStabilizer(stab1)
                stab2 = random_stabilizer(n)
                genstab2 = GeneralizedStabilizer(stab2)
                # apply some (repeated) non-Clifford operations to genstab1
                for _ in 1:rand(1:5)
                    i = rand(1:n)
                    nc = embed(n, i, pcT)
                    apply!(genstab1, nc)
                end
                # apply some (repeated) non-Clifford operations to genstab2
                for _ in 1:rand(1:5)
                    i = rand(1:n)
                    nc = embed(n, i, pcT)
                    apply!(genstab2, nc)
                end
                @test Operator(genstab1 ⊗ stab1) ≈ Operator(genstab1) ⊗ Operator(stab1)
                @test Operator(genstab2 ⊗ stab1) ≈ Operator(genstab2) ⊗ Operator(stab1)
                @test Operator(genstab1 ⊗ stab2) ≈ Operator(genstab1) ⊗ Operator(stab2)
                @test Operator(genstab2 ⊗ stab2) ≈ Operator(genstab2) ⊗ Operator(stab2)
                @test Operator(genstab1 ⊗ genstab2) ≈ Operator(genstab1) ⊗ Operator(genstab2)
                @test Operator(genstab2 ⊗ genstab1) ≈ Operator(genstab2) ⊗ Operator(genstab1)
                @test Operator(genstab1 ⊗ genstab1) ≈ Operator(genstab1) ⊗ Operator(genstab1)
                @test Operator(genstab2 ⊗ genstab2) ≈ Operator(genstab2) ⊗ Operator(genstab2)
            end
        end
    end

    @testset "Tensor products between paulichannels and paulis" begin
        num_trials = 3
        num_qubits = [2,3,4,5] # exclusively multi-qubit
        for n in num_qubits
            for repetition in 1:num_trials
                p = random_pauli(n)
                i = rand(1:n)
                nc = embed(n, i, pcT)
                @test Operator(nc ⊗ p) ≈ Operator(nc) ⊗ Operator(p)
            end
        end
    end

    @testset "smaller test redundant to the ones above" begin
        for n in 1:5
            for rep in 1:2
                s = random_stabilizer(n)
                g = GeneralizedStabilizer(s)
                apply!(g, embed(n, rand(1:n), pcT))
                p = random_pauli(n; realphase=true)
                gm, r = projectrand!(copy(g), p)

                rho = Operator(g)
                pqo = Operator(p)
                id = identityoperator(pqo)
                projp = (pqo+id)/2
                projm = (-pqo+id)/2

                @test projp+projm ≈ id

                rhom = projm*rho*projm'
                rhop = projp*rho*projp'

                #@test rhom + rhop ≈ rho

                @test (expect(p, g)+1)/2 ≈ tr(rhop)

                gm_notnorm, _ = QuantumClifford._projectrand_notnorm(copy(g), p, 0)
                @test (expect(p, g)+1)/2 ≈ tr(gm_notnorm)

                @test tr(rhop) ≈ tr(gm_notnorm)

                if r == 0x2
                    @test rhom / tr(rhom) ≈ Operator(gm)
                else
                    @test rhop / tr(rhop) ≈ Operator(gm)
                end
            end
        end
    end
end
