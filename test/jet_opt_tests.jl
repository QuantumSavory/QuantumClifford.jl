@testitem "JET optimization: projective measurements" tags=[:jet] begin
    using JET
    using QuantumClifford
    using Test

    state = one(MixedDestabilizer, 2)

    JET.@test_opt target_modules=(QuantumClifford,) projectY!(copy(state), 1; phases=true)
    JET.@test_opt target_modules=(QuantumClifford,) projectY!(copy(state), 1; phases=false)
    JET.@test_opt target_modules=(QuantumClifford,) projectZ!(copy(state), 1; phases=true)
    JET.@test_opt target_modules=(QuantumClifford,) projectZ!(copy(state), 1; phases=false)
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
