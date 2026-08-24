@testset "Stabilizer-state metrics match QuantumOpticsBase" begin
    using QuantumClifford
    import QuantumOpticsBase
    using Test

    function signed_stabilizer_group(state)
        group = [zero(PauliOperator, nqubits(state))]
        for generator in stabilizerview(state)
            products = [generator * element for element in group]
            append!(group, products)
        end
        group
    end

    function quantumoptics_density(state)
        density = exp2(-nqubits(state)) * sum(
            QuantumOpticsBase.Operator(pauli)
            for pauli in signed_stabilizer_group(state)
        )
        QuantumOpticsBase.dense(density)
    end

    function quantumoptics_expect(operator_state, state)
        operator = quantumoptics_density(operator_state)
        density = quantumoptics_density(state)
        QuantumOpticsBase.expect(operator, density)
    end

    function quantumoptics_expect(indices, operator_state, state)
        operator = quantumoptics_density(operator_state)
        density = quantumoptics_density(state)
        embedded = QuantumOpticsBase.embed(
            QuantumOpticsBase.basis(density),
            indices,
            operator,
        )
        QuantumOpticsBase.expect(embedded, density)
    end

    function quantumoptics_fidelity(state1, state2)
        density1 = quantumoptics_density(state1)
        density2 = quantumoptics_density(state2)
        QuantumOpticsBase.fidelity(density1, density2)
    end

    @testset "full-system expectation" begin
        pairs = (
            (MixedDestabilizer(S"Z"), MixedDestabilizer(S"Z")),
            (MixedDestabilizer(S"Z"), MixedDestabilizer(S"-Z")),
            (MixedDestabilizer(S"X"), MixedDestabilizer(S"-Y")),
            (MixedDestabilizer(bell()), MixedDestabilizer(S"Z_ _Z")),
            (MixedDestabilizer(S"-YY ZZ"), MixedDestabilizer(bell())),
            (maximally_mixed(1), maximally_mixed(1)),
            (maximally_mixed(2), maximally_mixed(2)),
            (MixedDestabilizer(S"Z_"), MixedDestabilizer(S"-Z_")),
            (MixedDestabilizer(S"Z_"), MixedDestabilizer(S"X_")),
        )

        for (operator_state, state) in pairs
            @test QuantumClifford.expect(operator_state, state) ≈
                real(quantumoptics_expect(operator_state, state)) atol=1e-12
        end
    end

    @testset "ordered subsystem expectation" begin
        asymmetric_operator = MixedDestabilizer(S"Z_ -_Z")
        asymmetric_state = MixedDestabilizer(S"Z__ _Z_ -__Z")
        cases = (
            (1, MixedDestabilizer(S"Z"), MixedDestabilizer(bell())),
            ([1, 3], MixedDestabilizer(bell()), MixedDestabilizer(ghz(3))),
            ([1, 3], asymmetric_operator, asymmetric_state),
            ([3, 1], asymmetric_operator, asymmetric_state),
            ([2], MixedDestabilizer(S"Z"), MixedDestabilizer(S"_Z_")),
        )

        for (indices, operator_state, state) in cases
            @test QuantumClifford.expect(indices, operator_state, state) ≈
                real(quantumoptics_expect(indices, operator_state, state)) atol=1e-12
        end
    end

    @testset "root fidelity" begin
        maximally_mixed_one = maximally_mixed(1)
        rank_one = MixedDestabilizer(S"Z_")
        pairs = (
            (MixedDestabilizer(S"Z"), MixedDestabilizer(S"Z")),
            (MixedDestabilizer(S"Z"), MixedDestabilizer(S"-Z")),
            (MixedDestabilizer(S"Z"), MixedDestabilizer(S"X")),
            (MixedDestabilizer(bell()), MixedDestabilizer(S"Z_ _Z")),
            (MixedDestabilizer(S"-YY ZZ"), MixedDestabilizer(bell())),
            (MixedDestabilizer(S"YY ZZ"), MixedDestabilizer(bell())),
            (MixedDestabilizer(S"Z"), maximally_mixed_one),
            (MixedDestabilizer(S"Z_ _Z"), rank_one),
            (MixedDestabilizer(S"X_ _X"), rank_one),
            (MixedDestabilizer(S"-Z_ _Z"), rank_one),
        )

        for (state1, state2) in pairs
            # QuantumOpticsBase evaluates this through dense matrix square roots.
            @test QuantumClifford.fidelity(state1, state2) ≈
                real(quantumoptics_fidelity(state1, state2)) atol=1e-7
        end
    end
end
