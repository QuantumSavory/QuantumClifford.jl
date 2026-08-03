@testitem "Stabilizer-state expectations" begin
    using LinearAlgebra
    using Random
    using StableRNGs
    using Test
    using QuantumClifford
    using QuantumInterface

    function signed_stabilizer_group(state)
        group = [zero(PauliOperator, nqubits(state))]
        for generator in stabilizerview(state)
            products = [generator * element for element in group]
            append!(group, products)
        end
        group
    end

    function expansion_expect(operator_state, state; indices=nothing)
        paulis = signed_stabilizer_group(operator_state)
        if !isnothing(indices)
            paulis = [embed(nqubits(state), indices, pauli) for pauli in paulis]
        end
        value = exp2(-nqubits(operator_state)) * sum(pauli -> expect(pauli, state), paulis)
        @test isreal(value)
        Float64(real(value))
    end

    mixed_state(rng, rank, qubits) = rank == 0 ? maximally_mixed(qubits) :
        MixedDestabilizer(random_stabilizer(rng, rank, qubits))

    @testset "pure states and signed generators" begin
        axes = (S"X", S"Y", S"Z")
        negative_axes = (S"-X", S"-Y", S"-Z")
        for (positive, negative) in zip(axes, negative_axes)
            @test expect(positive, positive) == 1.0
            @test expect(negative, negative) == 1.0
            @test expect(positive, negative) == 0.0
            @test expect(negative, positive) == 0.0
        end
        for first_axis in axes, second_axis in axes
            first_axis == second_axis && continue
            @test expect(first_axis, second_axis) ≈ 0.5
        end
        @test expect(S"-X", S"-Y") ≈ 0.5
        @test expect(S"X", S"-Y") ≈ 0.5

        bell_state = bell()
        product_state = S"Z_ _Z"
        @test expect(bell_state, product_state) ≈ 0.5
        @test expect(product_state, bell_state) ≈ 0.5
        @test expect(S"-YY ZZ", bell_state) == 1.0
        @test expect(S"YY ZZ", bell_state) == 0.0
    end

    @testset "rank-deficient states" begin
        maximally_mixed_one = maximally_mixed(1)
        maximally_mixed_two = maximally_mixed(2)
        rank_one = MixedDestabilizer(S"Z_")

        @test expect(S"Z", maximally_mixed_one) == 0.5
        @test expect(maximally_mixed_one, S"Z") == 0.5
        @test expect(maximally_mixed_one, maximally_mixed_one) == 0.5
        @test expect(maximally_mixed_two, maximally_mixed_two) == 0.25
        @test expect(rank_one, rank_one) == 0.5
        @test expect(MixedDestabilizer(S"Z_"), MixedDestabilizer(S"-Z_")) == 0.0
        @test expect(MixedDestabilizer(S"Z_"), MixedDestabilizer(S"X_")) == 0.25
    end

    @testset "redundant generators" begin
        redundant_pure = S"Z Z"
        pure = MixedDestabilizer(S"Z")
        @test expect(redundant_pure, pure) == 1.0
        @test expect(pure, redundant_pure) == 1.0
        @test expect(redundant_pure, redundant_pure) == 1.0
        @test fidelity(redundant_pure, pure) == 1.0

        redundant_mixed = S"Z_ Z_"
        mixed = MixedDestabilizer(S"Z_")
        pure_extension = S"Z_ _Z"
        @test expect(redundant_mixed, mixed) == 0.5
        @test expect(mixed, redundant_mixed) == 0.5
        @test expect(redundant_mixed, redundant_mixed) == 0.5
        @test fidelity(redundant_mixed, pure_extension) ≈ inv(sqrt(2))
        @test fidelity(pure_extension, redundant_mixed) ≈ inv(sqrt(2))
        @test_throws DomainError fidelity(redundant_mixed, mixed)

        redundant_bell = S"XX ZZ -YY"
        @test expect(redundant_bell, bell()) == 1.0
        @test expect(bell(), redundant_bell) == 1.0
        @test fidelity(redundant_bell, bell()) == 1.0
    end

    @testset "ordered subsystem embedding" begin
        ghz_state = ghz(3)
        @test expect(1, S"Z", bell()) == 0.5
        @test expect([1, 2], bell(), ghz_state) == 0.5
        @test expect((1, 3), bell(), ghz_state) == 0.5

        asymmetric_state = S"Z__ _Z_ -__Z"
        asymmetric_operator = S"Z_ -_Z"
        @test expect([1, 3], asymmetric_operator, asymmetric_state) == 1.0
        @test expect([3, 1], asymmetric_operator, asymmetric_state) == 0.0

        @test_throws DimensionMismatch expect([1], bell(), ghz_state)
        @test_throws ArgumentError expect([1, 1], bell(), ghz_state)
        @test_throws ArgumentError expect([0], S"Z", ghz_state)
        @test_throws ArgumentError expect([4], S"Z", ghz_state)
    end

    @testset "dimensions, representations, and immutability" begin
        @test_throws DimensionMismatch expect(S"Z", S"ZZ")

        pure = S"XX ZZ"
        pure_representations = (
            pure,
            Destabilizer(copy(pure)),
            MixedStabilizer(copy(pure)),
            MixedDestabilizer(copy(pure)),
            fastrow(MixedDestabilizer(copy(pure))),
            fastcolumn(MixedDestabilizer(copy(pure))),
        )
        for operator_state in pure_representations, state in pure_representations
            @test expect(operator_state, state) == 1.0
        end

        mixed = S"Z_"
        mixed_representations = (
            mixed,
            Destabilizer(copy(mixed)),
            MixedStabilizer(copy(mixed)),
            MixedDestabilizer(copy(mixed)),
            fastrow(MixedDestabilizer(copy(mixed))),
            fastcolumn(MixedDestabilizer(copy(mixed))),
        )
        for operator_state in mixed_representations, state in mixed_representations
            @test expect(operator_state, state) == 0.5
        end

        subsystem_state = canonicalize_noncomm(T"Z")
        subsystem_before = copy(subsystem_state)
        @test expect(subsystem_state, S"Z") == 1.0
        @test expect(S"Z", subsystem_state) == 1.0
        @test fidelity(subsystem_state, S"Z") == 1.0
        @test subsystem_state == subsystem_before

        operator_state = MixedDestabilizer(S"-YY ZZ")
        state = MixedDestabilizer(bell())
        operator_before = copy(operator_state)
        state_before = copy(state)
        @test expect(operator_state, state) == 1.0
        @test operator_state == operator_before
        @test state == state_before

        indexed_state = MixedDestabilizer(ghz(3))
        indexed_before = copy(indexed_state)
        @test expect([3, 1], S"X_ _Z", indexed_state) == 0.25
        @test indexed_state == indexed_before
    end

    @testset "independent signed-group expansion" begin
        rng = StableRNG(0x5a17)
        for qubits in 1:5, _ in 1:8
            operator_state = mixed_state(rng, rand(rng, 0:qubits), qubits)
            state = mixed_state(rng, rand(rng, 0:qubits), qubits)
            expected = expansion_expect(operator_state, state)
            value = expect(operator_state, state)
            reverse_value = expect(state, operator_state)
            @test value == reverse_value
            @test value ≈ expected
        end

        for state_qubits in 2:5, _ in 1:8
            operator_qubits = rand(rng, 1:state_qubits)
            operator_state = mixed_state(
                rng,
                rand(rng, 0:operator_qubits),
                operator_qubits,
            )
            state = mixed_state(rng, rand(rng, 0:state_qubits), state_qubits)
            indices = randperm(rng, state_qubits)[1:operator_qubits]
            @test expect(indices, operator_state, state) ==
                expansion_expect(operator_state, state; indices)
        end
    end

    @testset "root fidelity" begin
        @test QuantumClifford.fidelity === QuantumInterface.fidelity

        pure_pairs = (
            (S"Z", S"Z"),
            (S"Z", S"-Z"),
            (S"Z", S"X"),
            (bell(), S"Z_ _Z"),
            (bell(), S"-YY ZZ"),
        )
        for (state1, state2) in pure_pairs
            @test fidelity(state1, state2) == dot(state1, state2)
        end

        mixed = maximally_mixed(1)
        @test fidelity(S"Z", mixed) == fidelity(mixed, S"Z")
        @test fidelity(S"Z", mixed) ≈ inv(sqrt(2))
        @test abs2(fidelity(S"Z", mixed)) ≈ expect(S"Z", mixed)

        rank_one = MixedDestabilizer(S"Z_")
        @test fidelity(S"Z_ _Z", rank_one) ≈ inv(sqrt(2))
        @test abs2(fidelity(S"X_ _X", rank_one)) ==
            expect(S"X_ _X", rank_one) == 0.25
        @test fidelity(S"-Z_ _Z", rank_one) == 0.0

        @test_throws DomainError fidelity(maximally_mixed(1), maximally_mixed(1))
        @test_throws DomainError fidelity(
            MixedDestabilizer(S"Z_"),
            MixedDestabilizer(S"X_"),
        )
        @test_throws DimensionMismatch fidelity(S"Z", S"ZZ")

        pure_state = MixedDestabilizer(bell())
        mixed_state_value = MixedDestabilizer(S"ZZ")
        pure_before = copy(pure_state)
        mixed_before = copy(mixed_state_value)
        @test fidelity(pure_state, mixed_state_value) ≈ inv(sqrt(2))
        @test pure_state == pure_before
        @test mixed_state_value == mixed_before
    end
end
