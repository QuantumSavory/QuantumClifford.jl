@testitem "Word-boundary compaction after projectremoverand!" begin
    # Regression test for the word-boundary mismatch bug:
    # After projectremoverand! reduces a MixedDestabilizer across a UInt64
    # word boundary (e.g. 65→64 qubits), the internal xzs array must be
    # compacted so that `size(xzs,1) ÷ 2 == cld(nqubits, 64)`. Without
    # compaction, comm() misreads Z bits as X bits, causing expect() and
    # project!() to return silently wrong results.

    using QuantumClifford
    using QuantumClifford: MixedDestabilizer, projectremoverand!, tab,
        expect, comm, tensor, nqubits, stabilizerview

    # Sizes chosen to exercise word-boundary crossings:
    # 65→64 (2→1 words), 129→128 (3→2 words), and sizes that stay
    # within the same word count as controls.
    test_sizes = [1,2,10,63,64,65,127,128,129]

    @testset "xzs compaction invariant (n=$n)" for n in test_sizes
        n > 1 || continue # need at least 2 qubits for CZ + removal
        # Build n qubits in |+⟩ and apply CZ(1,2) to create a known stabilizer.
        plus = S"X"
        state = MixedDestabilizer(copy(plus))
        for i in 2:n
            state = tensor(state, MixedDestabilizer(copy(plus)))
        end
        apply!(state, sZCZ(1, 2))

        # Remove the last qubit via X-measurement.
        state, _ = projectremoverand!(state, projectX!, n)

        # Check the word-count invariant: physical words must match logical words.
        t = tab(state)
        physical_words = size(t.xzs, 1) ÷ 2
        logical_words  = cld(nqubits(state), 8 * sizeof(eltype(t.xzs)))
        @test physical_words == logical_words
    end

    @testset "expect correctness after word-boundary crossing (n=$n)" for n in [65, 129]
        # Build a path graph state P₃ on qubits 1-2-3, with extra qubits
        # in |+⟩ that push the total past a word boundary.
        plus = S"X"
        state = MixedDestabilizer(copy(plus))
        for i in 2:n
            state = tensor(state, MixedDestabilizer(copy(plus)))
        end
        apply!(state, sZCZ(1, 2))
        apply!(state, sZCZ(2, 3))

        # Remove qubits from the end until we cross the word boundary.
        for q in n:-1:(n == 65 ? 64 : 128)
            state, _ = projectremoverand!(state, projectX!, nqubits(state))
        end

        # X₁Z₂ is a stabilizer of the path graph P₃ — expect must be +1.
        nq = nqubits(state)
        op = PauliOperator(0x0,
            vcat([true, false], falses(nq - 2)),
            vcat([false, true], falses(nq - 2)))
        @test expect(op, state) == 1

        # Z₁X₂Z₃ is also a stabilizer of P₃.
        # X bits: qubit 2 only; Z bits: qubits 1 and 3 only (not 2, which would give Y₂).
        op2 = PauliOperator(0x0,
            vcat([false, true, false], falses(nq - 3)),
            vcat([true, false, true],  falses(nq - 3)))
        @test expect(op2, state) == 1
    end

    @testset "comm consistency after word-boundary crossing" begin
        # Build 65-qubit state, apply CZ(1,2), remove qubit 65 → 64 qubits.
        plus = S"X"
        state = MixedDestabilizer(copy(plus))
        for i in 2:65
            state = tensor(state, MixedDestabilizer(copy(plus)))
        end
        apply!(state, sZCZ(1, 2))
        state, _ = projectremoverand!(state, projectX!, 65)

        # Build the same 64-qubit state directly (no word-boundary crossing).
        state_direct = MixedDestabilizer(copy(plus))
        for i in 2:64
            state_direct = tensor(state_direct, MixedDestabilizer(copy(plus)))
        end
        apply!(state_direct, sZCZ(1, 2))

        # X₁Z₂ must commute with all stabilizer rows in both states.
        nq = nqubits(state)
        op = PauliOperator(0x0,
            vcat([true, false], falses(nq - 2)),
            vcat([false, true], falses(nq - 2)))

        stab_crossed = stabilizerview(state)
        stab_direct  = stabilizerview(state_direct)
        for i in 1:length(stab_crossed)
            @test comm(op, stab_crossed, i) == 0x00
        end
        for i in 1:length(stab_direct)
            @test comm(op, stab_direct, i) == 0x00
        end
    end

    @testset "repeated removal across multiple word boundaries" begin
        # Start at 130 qubits (3 words), remove down to 63 (1 word).
        # This crosses two word boundaries: 3→2 words and 2→1 words.
        plus = S"X"
        state = MixedDestabilizer(copy(plus))
        for i in 2:130
            state = tensor(state, MixedDestabilizer(copy(plus)))
        end
        apply!(state, sZCZ(1, 2))

        for q in 130:-1:64
            state, _ = projectremoverand!(state, projectX!, nqubits(state))
        end

        # After all removals, the state should have 63 qubits and 1 word.
        @test nqubits(state) == 63
        t = tab(state)
        @test size(t.xzs, 1) ÷ 2 == 1

        # X₁Z₂ is still a stabilizer (CZ on qubits 1-2 was never removed).
        nq = nqubits(state)
        op = PauliOperator(0x0,
            vcat([true, false], falses(nq - 2)),
            vcat([false, true], falses(nq - 2)))
        @test expect(op, state) == 1
    end
end
