struct ZeroBasedIndices <: AbstractVector{Int}
    data::Vector{Int}
end
Base.size(indices::ZeroBasedIndices) = size(indices.data)
Base.axes(indices::ZeroBasedIndices) = (0:(length(indices.data) - 1),)
Base.getindex(indices::ZeroBasedIndices, index::Int) = indices.data[index + 1]

@testset "Stabilizer-state expectations" begin
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

    function repack_state(state, ::Type{T}) where {T<:Unsigned}
        source_state = MixedDestabilizer(state)
        source = QuantumClifford.tab(source_state)
        target = zero(
            QuantumClifford.Tableau{Vector{UInt8},Matrix{T}},
            length(source),
            nqubits(source),
        )
        phases(target) .= UInt8.(phases(source))
        for row in eachindex(source), qubit in 1:nqubits(source)
            target[row, qubit] = source[row, qubit]
        end
        MixedDestabilizer(target, rank(source_state))
    end

    @testset "pure states and signed generators" begin
        axes = (
            MixedDestabilizer(S"X"),
            MixedDestabilizer(S"Y"),
            MixedDestabilizer(S"Z"),
        )
        negative_axes = (
            MixedDestabilizer(S"-X"),
            MixedDestabilizer(S"-Y"),
            MixedDestabilizer(S"-Z"),
        )
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
        @test expect(MixedDestabilizer(S"-X"), MixedDestabilizer(S"-Y")) ≈ 0.5
        @test expect(MixedDestabilizer(S"X"), MixedDestabilizer(S"-Y")) ≈ 0.5

        bell_state = MixedDestabilizer(bell())
        product_state = MixedDestabilizer(S"Z_ _Z")
        @test expect(bell_state, product_state) ≈ 0.5
        @test expect(product_state, bell_state) ≈ 0.5
        @test expect(MixedDestabilizer(S"-YY ZZ"), bell_state) == 1.0
        @test expect(MixedDestabilizer(S"YY ZZ"), bell_state) == 0.0
    end

    @testset "rank-deficient states" begin
        maximally_mixed_one = maximally_mixed(1)
        maximally_mixed_two = maximally_mixed(2)
        rank_one = MixedDestabilizer(S"Z_")

        pure_z = MixedDestabilizer(S"Z")
        @test expect(pure_z, maximally_mixed_one) == 0.5
        @test expect(maximally_mixed_one, pure_z) == 0.5
        @test expect(maximally_mixed_one, maximally_mixed_one) == 0.5
        @test expect(maximally_mixed_two, maximally_mixed_two) == 0.25
        @test expect(rank_one, rank_one) == 0.5
        @test expect(MixedDestabilizer(S"Z_"), MixedDestabilizer(S"-Z_")) == 0.0
        @test expect(MixedDestabilizer(S"Z_"), MixedDestabilizer(S"X_")) == 0.25
    end

    @testset "redundant generators" begin
        redundant_pure = MixedDestabilizer(S"Z Z")
        pure = MixedDestabilizer(S"Z")
        @test expect(redundant_pure, pure) == 1.0
        @test expect(pure, redundant_pure) == 1.0
        @test expect(redundant_pure, redundant_pure) == 1.0
        @test fidelity(redundant_pure, pure) == 1.0

        redundant_mixed = MixedDestabilizer(S"Z_ Z_")
        mixed = MixedDestabilizer(S"Z_")
        pure_extension = MixedDestabilizer(S"Z_ _Z")
        @test expect(redundant_mixed, mixed) == 0.5
        @test expect(mixed, redundant_mixed) == 0.5
        @test expect(redundant_mixed, redundant_mixed) == 0.5
        @test fidelity(redundant_mixed, pure_extension) ≈ inv(sqrt(2))
        @test fidelity(pure_extension, redundant_mixed) ≈ inv(sqrt(2))
        @test_throws DomainError fidelity(redundant_mixed, mixed)

        redundant_bell = MixedDestabilizer(S"XX ZZ -YY")
        bell_state = MixedDestabilizer(bell())
        @test expect(redundant_bell, bell_state) == 1.0
        @test expect(bell_state, redundant_bell) == 1.0
        @test fidelity(redundant_bell, bell_state) == 1.0
    end

    @testset "ordered subsystem embedding" begin
        ghz_state = MixedDestabilizer(ghz(3))
        @test expect(
            1,
            MixedDestabilizer(S"Z"),
            MixedDestabilizer(bell()),
        ) == 0.5
        @test expect([1, 2], MixedDestabilizer(bell()), ghz_state) == 0.5
        @test expect((1, 3), MixedDestabilizer(bell()), ghz_state) == 0.5

        asymmetric_state = MixedDestabilizer(S"Z__ _Z_ -__Z")
        asymmetric_operator = MixedDestabilizer(S"Z_ -_Z")
        @test expect([1, 3], asymmetric_operator, asymmetric_state) == 1.0
        @test expect([3, 1], asymmetric_operator, asymmetric_state) == 0.0
        zero_based_indices = ZeroBasedIndices([1, 3])
        @test_throws ArgumentError expect(
            zero_based_indices,
            asymmetric_operator,
            asymmetric_state,
        )
        @test zero_based_indices.data == [1, 3]

        bell_state = MixedDestabilizer(bell())
        z_state = MixedDestabilizer(S"Z")
        @test_throws DimensionMismatch expect([1], bell_state, ghz_state)
        @test_throws ArgumentError expect([1, 1], bell_state, ghz_state)
        @test_throws ArgumentError expect([0], z_state, ghz_state)
        @test_throws ArgumentError expect([4], z_state, ghz_state)
        @test_throws DimensionMismatch expect(1, bell_state, ghz_state)
        @test_throws ArgumentError expect(0, z_state, ghz_state)
        @test_throws ArgumentError expect(4, z_state, ghz_state)
    end

    @testset "dimensions, representations, and immutability" begin
        @test_throws DimensionMismatch expect(
            MixedDestabilizer(S"Z"),
            MixedDestabilizer(S"ZZ"),
        )

        pure = S"XX ZZ"
        pure_representations = (
            MixedDestabilizer(copy(pure)),
            fastrow(MixedDestabilizer(copy(pure))),
            fastcolumn(MixedDestabilizer(copy(pure))),
        )
        for operator_state in pure_representations, state in pure_representations
            @test expect(operator_state, state) == 1.0
        end

        mixed = S"Z_"
        mixed_representations = (
            MixedDestabilizer(copy(mixed)),
            fastrow(MixedDestabilizer(copy(mixed))),
            fastcolumn(MixedDestabilizer(copy(mixed))),
        )
        for operator_state in mixed_representations, state in mixed_representations
            @test expect(operator_state, state) == 0.5
        end

        subsystem_state = canonicalize_noncomm(T"Z")
        subsystem_before = copy(subsystem_state)
        converted_subsystem_state = MixedDestabilizer(stabilizerview(subsystem_state))
        pure_z = MixedDestabilizer(S"Z")
        @test expect(converted_subsystem_state, pure_z) == 1.0
        @test expect(pure_z, converted_subsystem_state) == 1.0
        @test fidelity(converted_subsystem_state, pure_z) == 1.0
        @test subsystem_state == subsystem_before

        operator_state = MixedDestabilizer(S"-YY ZZ")
        state = MixedDestabilizer(bell())
        operator_before = copy(operator_state)
        state_before = copy(state)
        @test expect(operator_state, state) == 1.0
        @test operator_state == operator_before
        @test state == state_before

        indexed_operator = MixedDestabilizer(S"X_ _Z")
        indexed_state = MixedDestabilizer(ghz(3))
        indexed_operator_before = copy(indexed_operator)
        indexed_before = copy(indexed_state)
        @test expect(
            [3, 1],
            indexed_operator,
            indexed_state,
        ) == 0.25
        @test indexed_operator == indexed_operator_before
        @test indexed_state == indexed_before
    end

    @testset "packed word representations" begin
        for word_type in (UInt8, UInt16, UInt32, UInt64)
            packed_z = repack_state(MixedDestabilizer(S"Z"), word_type)
            packed_maximally_mixed = repack_state(maximally_mixed(1), word_type)
            @test expect(packed_z, packed_maximally_mixed) == 0.5
            @test expect((1,), packed_z, packed_maximally_mixed) == 0.5
            @test fidelity(packed_z, packed_maximally_mixed) ≈ inv(sqrt(2))

            packed_mixed = repack_state(MixedDestabilizer(S"Z_"), word_type)
            @test expect(packed_mixed, packed_mixed) == 0.5

            packed_pure = repack_state(MixedDestabilizer(bell()), word_type)
            @test expect(packed_pure, packed_pure) == 1.0
            @test fidelity(packed_pure, packed_pure) == 1.0
        end

        long_mixed = MixedDestabilizer([embed(65, 65, P"Z")])
        for word_type in (UInt8, UInt16, UInt32, UInt64)
            packed_long_mixed = repack_state(long_mixed, word_type)
            packed_z = repack_state(MixedDestabilizer(S"Z"), word_type)
            @test expect(packed_long_mixed, packed_long_mixed) == exp2(-64)
            @test expect((65,), packed_z, packed_long_mixed) == 1.0
            @test expect((1,), packed_z, packed_long_mixed) == 0.5
        end
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
            (MixedDestabilizer(S"Z"), MixedDestabilizer(S"Z")),
            (MixedDestabilizer(S"Z"), MixedDestabilizer(S"-Z")),
            (MixedDestabilizer(S"Z"), MixedDestabilizer(S"X")),
            (MixedDestabilizer(bell()), MixedDestabilizer(S"Z_ _Z")),
            (MixedDestabilizer(bell()), MixedDestabilizer(S"-YY ZZ")),
        )
        for (state1, state2) in pure_pairs
            @test fidelity(state1, state2) == dot(state1, state2)
        end

        mixed = maximally_mixed(1)
        pure_z = MixedDestabilizer(S"Z")
        @test fidelity(pure_z, mixed) == fidelity(mixed, pure_z)
        @test fidelity(pure_z, mixed) ≈ inv(sqrt(2))
        @test abs2(fidelity(pure_z, mixed)) ≈ expect(pure_z, mixed)

        rank_one = MixedDestabilizer(S"Z_")
        pure_extension = MixedDestabilizer(S"Z_ _Z")
        mutually_unbiased = MixedDestabilizer(S"X_ _X")
        @test fidelity(pure_extension, rank_one) ≈ inv(sqrt(2))
        @test abs2(fidelity(mutually_unbiased, rank_one)) ==
            expect(mutually_unbiased, rank_one) == 0.25
        @test fidelity(MixedDestabilizer(S"-Z_ _Z"), rank_one) == 0.0

        @test_throws DomainError fidelity(maximally_mixed(1), maximally_mixed(1))
        @test_throws DomainError fidelity(
            MixedDestabilizer(S"Z_"),
            MixedDestabilizer(S"X_"),
        )
        @test_throws DimensionMismatch fidelity(
            MixedDestabilizer(S"Z"),
            MixedDestabilizer(S"ZZ"),
        )

        pure_state = MixedDestabilizer(bell())
        mixed_state_value = MixedDestabilizer(S"ZZ")
        pure_before = copy(pure_state)
        mixed_before = copy(mixed_state_value)
        @test fidelity(pure_state, mixed_state_value) ≈ inv(sqrt(2))
        @test pure_state == pure_before
        @test mixed_state_value == mixed_before
    end
end
