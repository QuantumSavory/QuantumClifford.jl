@testitem "Word-boundary compaction after qubit removal" begin
    using QuantumClifford
    using QuantumClifford: MixedDestabilizer, projectremoverand!, remove_column!,
        tab, expect, comm, tensor, nqubits, stabilizerview

    function plus_destabilizer(n)
        plus = S"X"
        state = MixedDestabilizer(copy(plus))
        for _ in 2:n
            state = tensor(state, MixedDestabilizer(copy(plus)))
        end
        state
    end

    function pauli_on(n; xs=(), zs=())
        xbits = falses(n)
        zbits = falses(n)
        for i in xs
            xbits[i] = true
        end
        for i in zs
            zbits[i] = true
        end
        PauliOperator(0x0, xbits, zbits)
    end

    function has_compact_xzs(state)
        t = tab(state)
        size(t.xzs, 1) ÷ 2 == cld(nqubits(state), 8 * sizeof(eltype(t.xzs)))
    end

    invariant_sizes = [1, 2, 10, 63, 64, 65, 127, 128, 129]

    @testset "xzs word-count invariant (n=$n)" for n in invariant_sizes
        state = plus_destabilizer(n)
        n > 1 && apply!(state, sZCZ(1, 2))

        state, _ = projectremoverand!(state, projectX!, n)
        @test has_compact_xzs(state)
    end

    @testset "remove_column! compacts direct callers" begin
        state = plus_destabilizer(65)
        removed = remove_column!(copy(stabilizerview(state)), 65)

        @test nqubits(removed) == 64
        @test has_compact_xzs(removed)
        @test comm(pauli_on(64; xs=(1,)), removed, 1) in (0x00, 0x01)
    end

    @testset "expect correctness after word-boundary crossing (n=$n)" for n in [65, 129]
        state = plus_destabilizer(n)
        apply!(state, sZCZ(1, 2))
        apply!(state, sZCZ(2, 3))

        for _ in 1:2
            state, _ = projectremoverand!(state, projectX!, nqubits(state))
        end

        nq = nqubits(state)
        @test expect(pauli_on(nq; xs=(1,), zs=(2,)), state) == 1
        @test expect(pauli_on(nq; xs=(2,), zs=(1, 3)), state) == 1
    end

    @testset "comm consistency after word-boundary crossing" begin
        crossed = plus_destabilizer(65)
        apply!(crossed, sZCZ(1, 2))
        crossed, _ = projectremoverand!(crossed, projectX!, 65)

        direct = plus_destabilizer(64)
        apply!(direct, sZCZ(1, 2))

        op = pauli_on(64; xs=(1,), zs=(2,))
        @test comm(op, stabilizerview(crossed)) ==
            comm(op, stabilizerview(direct)) ==
            zeros(UInt8, length(stabilizerview(crossed)))
    end

    @testset "non-tail removal across word boundary" begin
        state = plus_destabilizer(65)
        apply!(state, sZCZ(64, 65))

        state, _ = projectremoverand!(state, projectX!, 1)

        @test nqubits(state) == 64
        @test has_compact_xzs(state)
        @test expect(pauli_on(64; xs=(63,), zs=(64,)), state) == 1
        @test expect(pauli_on(64; xs=(64,), zs=(63,)), state) == 1
    end

    @testset "repeated removal across multiple word boundaries" begin
        state = plus_destabilizer(130)
        apply!(state, sZCZ(1, 2))

        while nqubits(state) > 63
            state, _ = projectremoverand!(state, projectX!, nqubits(state))
        end

        @test nqubits(state) == 63
        @test has_compact_xzs(state)
        @test expect(pauli_on(63; xs=(1,), zs=(2,)), state) == 1
    end
end
