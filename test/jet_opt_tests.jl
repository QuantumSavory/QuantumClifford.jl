@testitem "JET optimization: projective measurements" tags=[:jet] begin
    using JET
    using QuantumClifford
    using Test

    state = one(MixedDestabilizer, 2)

    JET.@test_opt target_modules=(QuantumClifford,) projectX!(copy(state), 1; phases=true)
    JET.@test_opt target_modules=(QuantumClifford,) projectX!(copy(state), 1; phases=false)
    JET.@test_opt target_modules=(QuantumClifford,) projectY!(copy(state), 1; phases=true)
    JET.@test_opt target_modules=(QuantumClifford,) projectY!(copy(state), 1; phases=false)
    JET.@test_opt target_modules=(QuantumClifford,) projectZ!(copy(state), 1; phases=true)
    JET.@test_opt target_modules=(QuantumClifford,) projectZ!(copy(state), 1; phases=false)
end

@testitem "JET optimization: Pauli primitives and constructors" tags=[:jet] begin
    using JET
    using QuantumClifford
    using Test

    pauli = P"XZ"
    other_pauli = P"ZY"

    JET.@test_opt target_modules=(QuantumClifford,) comm(pauli, other_pauli)
    JET.@test_opt target_modules=(QuantumClifford,) pauli * other_pauli
    JET.@test_opt target_modules=(QuantumClifford,) pauli ⊗ other_pauli
    JET.@test_opt target_modules=(QuantumClifford,) QuantumClifford.mul_left!(copy(pauli), other_pauli)
    JET.@test_opt target_modules=(QuantumClifford,) random_pauli(4)
    JET.@test_opt target_modules=(QuantumClifford,) random_stabilizer(4)
    JET.@test_opt target_modules=(QuantumClifford,) random_destabilizer(4)
    JET.@test_opt target_modules=(QuantumClifford,) random_clifford(4)
end

@testitem "JET optimization: Clifford and state application" tags=[:jet] begin
    using JET
    using QuantumClifford
    using Test

    bell_state = bell()
    state = S"-XXX
              +ZZI
              -IZZ"
    mixed_bell = MixedDestabilizer(bell_state)

    JET.@test_opt target_modules=(QuantumClifford,) tensor(bell_state, bell_state)
    JET.@test_opt target_modules=(QuantumClifford,) tCNOT * bell_state
    JET.@test_opt target_modules=(QuantumClifford,) apply!(copy(bell_state), tCNOT)
    JET.@test_opt target_modules=(QuantumClifford,) apply!(copy(state), [1, 2], tCNOT)
    JET.@test_opt target_modules=(QuantumClifford,) apply!(copy(mixed_bell), sHadamard(1))
    JET.@test_opt target_modules=(QuantumClifford,) apply!(copy(mixed_bell), sCNOT(1, 2))
end

@testitem "JET optimization: tableau lifecycle" tags=[:jet] begin
    using JET
    using QuantumClifford
    using Test

    state = S"-XXX
              +ZZI
              -IZZ"
    canonical_state = canonicalize!(copy(state))
    mixed_state = MixedDestabilizer(state)

    JET.@test_opt target_modules=(QuantumClifford,) canonicalize!(copy(state))
    JET.@test_opt target_modules=(QuantumClifford,) generate!(P"XYY", canonical_state)
    JET.@test_opt target_modules=(QuantumClifford,) project!(copy(mixed_state), P"ZII")
    JET.@test_opt target_modules=(QuantumClifford,) expect(P"ZZ_", mixed_state)
    JET.@test_opt target_modules=(QuantumClifford,) traceout!(copy(mixed_state), [1])
    JET.@test_opt target_modules=(QuantumClifford,) ptrace(mixed_state, [1])
end

@testitem "JET optimization: qubit reset" tags=[:jet] begin
    using JET
    using QuantumClifford
    using Test

    state = MixedDestabilizer(bell())

    JET.@test_opt target_modules=(QuantumClifford,) reset_qubits!(
        copy(state), S"Z", [1]
    )
end

@testitem "JET optimization: Monte Carlo trajectories" tags=[:jet] begin
    using JET
    using QuantumClifford
    using Test

    initial_state = MixedDestabilizer(bell())
    circuit = [sHadamard(1), sCNOT(1, 2), sMZ(1)]

    # The heterogeneous circuit dispatches through abstract `applywstatus!` calls.
    JET.@test_opt target_modules=(QuantumClifford,) broken=true mctrajectories(
        initial_state, circuit; trajectories=2
    )
end

@testitem "JET optimization: naive encoding circuit" tags=[:jet] begin
    using JET
    using QuantumClifford
    using QuantumClifford.ECC: Steane7, naive_encoding_circuit
    using Test

    JET.@test_opt target_modules=(QuantumClifford, QuantumClifford.ECC) naive_encoding_circuit(Steane7())
end

@testitem "JET optimization: table decoders" tags=[:jet] begin
    using JET
    using QuantumClifford
    using QuantumClifford.ECC: CSSTableDecoder, CommutationCheckECCSetup, Shor9, Steane7, TableDecoder, evaluate_decoder
    using Test

    table_decoder = TableDecoder(Shor9())
    css_table_decoder = CSSTableDecoder(Steane7())

    @test isconcretetype(typeof(table_decoder))
    @test isconcretetype(typeof(css_table_decoder))
    @test all(isconcretetype, fieldtypes(typeof(table_decoder)))
    @test all(isconcretetype, fieldtypes(typeof(css_table_decoder)))

    table_from_fields = TableDecoder(
        table_decoder.H,
        table_decoder.faults_matrix,
        Int32(table_decoder.n),
        Int32(table_decoder.s),
        Int32(table_decoder.k),
        Int32(table_decoder.error_weight),
        table_decoder.lookup_table,
    )
    css_from_fields = CSSTableDecoder(
        css_table_decoder.H,
        css_table_decoder.faults_matrix,
        Int32(css_table_decoder.n),
        Int32(css_table_decoder.s),
        Int32(css_table_decoder.k),
        Int32(css_table_decoder.cx),
        Int32(css_table_decoder.cz),
        Int32(css_table_decoder.error_weight),
        css_table_decoder.tabledecoderx,
        css_table_decoder.tabledecoderz,
    )
    @test typeof(table_from_fields) === typeof(table_decoder)
    @test typeof(css_from_fields) === typeof(css_table_decoder)

    setup = CommutationCheckECCSetup(0.001)
    JET.@test_opt target_modules=(QuantumClifford, QuantumClifford.ECC) evaluate_decoder(table_decoder, setup, 32)
    JET.@test_opt target_modules=(QuantumClifford, QuantumClifford.ECC) evaluate_decoder(css_table_decoder, setup, 32)
end
