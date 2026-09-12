@testitem "decoder pipeline: decode!" begin
    using QuantumClifford
    using QuantumClifford: falses, trues

    @eval import QuantumClifford: ECC
    const E = ECC

    c = E.Steane7()
    s = E.code_s(c)
    n = E.code_n(c)

    td = E.TableDecoder(c)
    csstd = E.CSSTableDecoder(c)

    # decode! success case for TableDecoder
    syndrome = falses(s)
    expected = E.decode(td, syndrome)
    out = trues(2n)
    ret = E.decode!(out, td, syndrome)
    @test ret === out
    @test out == expected

    # decode! success case for CSSTableDecoder
    out2 = trues(2n)
    ret2 = E.decode!(out2, csstd, syndrome)
    @test ret2 === out2
    @test out2 == E.decode(csstd, syndrome)

    # decode! where syndrome is not in the lookup table (returns nothing)
    missing_syndrome = Bool[0, 1, 0, 1, 0, 0]
    @test isnothing(E.decode!(trues(2n), td, missing_syndrome))

    # ClassicalTableDecoder with custom matrix
    custom_H = Bool[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 1]
    ctd = E.ClassicalTableDecoder(custom_H)
    Hrows, Hcols = size(custom_H)
    # [1,0,0,0,1] is in the lookup table (column 1)
    @test !isnothing(E.decode!(trues(Hcols), ctd, Bool[1,0,0,0,1]))
    # [0,0,0,0,1] is missing from lookup table
    @test isnothing(E.decode!(trues(Hcols), ctd, Bool[0,0,0,0,1]))

    # Wrong-dimension errors for TableDecoder
    @test_throws ArgumentError E.decode!(trues(0), td, falses(s))

    # Wrong-dimension errors for ClassicalTableDecoder
    @test_throws ArgumentError E.decode!(trues(0), ctd, falses(Hrows))

    # Wrong-dimension errors for CSSTableDecoder
    @test_throws ArgumentError E.decode!(trues(0), csstd, falses(s))
    @test_throws ArgumentError E.decode!(trues(n), csstd, trues(s+1))
end

@testitem "decoder pipeline: batchdecode and batchdecode!" begin
    using QuantumClifford
    using QuantumClifford: falses, trues

    @eval import QuantumClifford: ECC
    const E = ECC

    c = E.Steane7()
    s = E.code_s(c)
    n = E.code_n(c)

    td = E.TableDecoder(c)
    csstd = E.CSSTableDecoder(c)
    Hx = E.parity_matrix_x(c)
    cx = size(Hx, 1)
    ctd = E.ClassicalTableDecoder(Matrix{Bool}(Hx))

    syndromes = falses(3, s)
    syndromes[2, 1] = true
    syndromes[3, 2] = true

    # batchdecode generic
    expected_batch = E.batchdecode(td, syndromes)
    @test size(expected_batch) == (3, 2n)

    # batchdecode! for TableDecoder (serial)
    out_batch = trues(3, 2n)
    ret_batch = E.batchdecode!(out_batch, td, syndromes)
    @test ret_batch === out_batch
    @test out_batch == expected_batch

    # batchdecode! for CSSTableDecoder (serial)
    out_css = trues(3, 2n)
    ret_css = E.batchdecode!(out_css, csstd, syndromes)
    @test ret_css === out_css
    @test out_css == E.batchdecode(csstd, syndromes)

    # ClassicalTableDecoder batchdecode! (serial)
    synd_x = syndromes[:, 1:cx]
    expected_ctd = vcat(
        reshape(E.decode(ctd, synd_x[1,:]), 1, n),
        reshape(E.decode(ctd, synd_x[2,:]), 1, n),
        reshape(E.decode(ctd, synd_x[3,:]), 1, n),
    )
    out_ctd = trues(3, n)
    ret_ctd = E.batchdecode!(out_ctd, ctd, synd_x)
    @test ret_ctd === out_ctd
    @test out_ctd == expected_ctd

    # decode_batch alias
    @test E.decode_batch(td, syndromes) == expected_batch

    # decode_batch! alias
    @test E.decode_batch!(trues(3, 2n), td, syndromes) == expected_batch
    @test E.decode_batch!(Matrix{Bool}(undef,3,2n), td, syndromes; threaded=true) == expected_batch

    # Wrong-dimension errors for TableDecoder batchdecode!
    @test_throws ArgumentError E.batchdecode!(trues(2,2n), td, trues(3,s+1))
    @test_throws ArgumentError E.batchdecode!(trues(1,2n), td, syndromes)
    @test_throws ArgumentError E.batchdecode!(trues(2,1), td, syndromes)

    # Wrong-dimension errors for ClassicalTableDecoder batchdecode!
    @test_throws ArgumentError E.batchdecode!(trues(2,n), ctd, trues(3,cx+1))
    @test_throws ArgumentError E.batchdecode!(trues(1,n), ctd, synd_x)
    @test_throws ArgumentError E.batchdecode!(trues(2,1), ctd, synd_x)

    # Wrong-dimension errors for CSSTableDecoder batchdecode!
    @test_throws ArgumentError E.batchdecode!(trues(2,2n), csstd, trues(3,s+1))
    @test_throws ArgumentError E.batchdecode!(trues(1,2n), csstd, syndromes)
    @test_throws ArgumentError E.batchdecode!(trues(2,1), csstd, syndromes)

    # Wrong-dimension errors for generic batchdecode
    @test_throws ArgumentError E.batchdecode(td, trues(2,s+1))
end

@testitem "decoder pipeline: batchdecode! threaded" begin
    using QuantumClifford
    using QuantumClifford: falses, trues

    @eval import QuantumClifford: ECC
    const E = ECC

    c = E.Steane7()
    s = E.code_s(c)
    n = E.code_n(c)

    td = E.TableDecoder(c)
    csstd = E.CSSTableDecoder(c)
    Hx = E.parity_matrix_x(c)
    cx = size(Hx, 1)
    ctd = E.ClassicalTableDecoder(Matrix{Bool}(Hx))

    syndromes = falses(3, s)
    syndromes[2, 1] = true
    syndromes[3, 2] = true

    expected_batch = E.batchdecode(td, syndromes)

    # Threaded with BitMatrix (should warn and fall back to serial)
    @test_logs (:warn,) E.batchdecode!(falses(3,2n), td, syndromes; threaded=true)

    # Threaded with Matrix{Bool}
    out_threaded = Matrix{Bool}(undef, 3, 2n)
    ret_threaded = E.batchdecode!(out_threaded, td, syndromes; threaded=true)
    @test ret_threaded === out_threaded
    @test out_threaded == expected_batch

    # CSSTableDecoder threaded with BitMatrix
    @test_logs (:warn,) E.batchdecode!(falses(3,2n), csstd, syndromes; threaded=true)

    # CSSTableDecoder threaded with Matrix{Bool}
    out_css_threaded = Matrix{Bool}(undef, 3, 2n)
    ret_css_threaded = E.batchdecode!(out_css_threaded, csstd, syndromes; threaded=true)
    @test ret_css_threaded === out_css_threaded
    @test out_css_threaded == E.batchdecode(csstd, syndromes)

    # ClassicalTableDecoder threaded with BitMatrix
    synd_x = syndromes[:, 1:cx]
    expected_ctd = vcat(
        reshape(E.decode(ctd, synd_x[1,:]), 1, n),
        reshape(E.decode(ctd, synd_x[2,:]), 1, n),
        reshape(E.decode(ctd, synd_x[3,:]), 1, n),
    )
    @test_logs (:warn,) E.batchdecode!(falses(3,n), ctd, synd_x; threaded=true)

    # ClassicalTableDecoder threaded with Matrix{Bool}
    out_ctd_threaded = Matrix{Bool}(undef, 3, n)
    ret_ctd_threaded = E.batchdecode!(out_ctd_threaded, ctd, synd_x; threaded=true)
    @test ret_ctd_threaded === out_ctd_threaded
    @test out_ctd_threaded == expected_ctd
end

@testitem "decoder pipeline: generic AbstractSyndromeDecoder fallback" begin
    using QuantumClifford
    using QuantumClifford: falses, trues

    @eval import QuantumClifford: ECC
    const E = ECC

    struct GenericFallbackDecoder <: E.AbstractSyndromeDecoder
        H
        n::Int
    end
    E.parity_checks(d::GenericFallbackDecoder) = d.H
    function E.decode(d::GenericFallbackDecoder, syndrome_sample)
        return falses(2d.n)
    end

    c = E.Steane7()
    H = E.parity_checks(c)
    s = E.code_s(c)
    n = E.code_n(c)
    d = GenericFallbackDecoder(H, n)

    # generic decode! (success case)
    out = trues(2n)
    ret = E.decode!(out, d, falses(s))
    @test ret === out
    @test out == falses(2n)

    # generic batchdecode
    syndromes = falses(3, s)
    syndromes[2, 1] = true
    result = E.batchdecode(d, syndromes)
    @test result == falses(3, 2n)

    # generic batchdecode! serial
    out_batch = trues(3, 2n)
    ret_batch = E.batchdecode!(out_batch, d, syndromes)
    @test ret_batch === out_batch
    @test out_batch == falses(3, 2n)

    # generic batchdecode! threaded (no specialized batchdecode! -> generic fallback with warning)
    @test_logs (:warn,) E.batchdecode!(Matrix{Bool}(undef,3,2n), d, syndromes; threaded=true)

    # Wrong-dimension errors for generic batchdecode
    @test_throws ArgumentError E.batchdecode(d, trues(2,s+1))

    # Wrong-dimension errors for generic batchdecode!
    @test_throws ArgumentError E.batchdecode!(trues(3,2n), d, trues(2,s+1))
    @test_throws ArgumentError E.batchdecode!(trues(2,2n), d, syndromes)
    @test_throws ArgumentError E.batchdecode!(trues(3,1), d, syndromes)
end

@testitem "decoder pipeline: decode_batch aliases" begin
    using QuantumClifford
    using QuantumClifford: falses, trues

    @eval import QuantumClifford: ECC
    const E = ECC

    c = E.Steane7()
    s = E.code_s(c)
    n = E.code_n(c)
    td = E.TableDecoder(c)

    syndromes = falses(3, s)
    syndromes[2, 1] = true

    expected = E.batchdecode(td, syndromes)
    @test E.decode_batch(td, syndromes) == expected
    @test E.decode_batch!(trues(3, 2n), td, syndromes) == expected
    @test E.decode_batch!(Matrix{Bool}(undef,3,2n), td, syndromes; threaded=true) == expected
end

@testitem "decoder pipeline: ClassicalTableDecoder edge cases" begin
    using QuantumClifford
    using QuantumClifford: falses, trues

    @eval import QuantumClifford: ECC
    const E = ECC

    custom_H = Bool[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 1]
    ctd = E.ClassicalTableDecoder(custom_H)
    Hrows, Hcols = size(custom_H)

    # Zero syndrome should return zero correction
    @test E.decode!(trues(Hcols), ctd, falses(Hrows)) == falses(Hcols)

    # Syndrome [1,0,0,0,1] is in table (column 1) -> correct qubit 1
    result1 = E.decode!(trues(Hcols), ctd, Bool[1,0,0,0,1])
    @test result1[1] == 1
    @test sum(result1) == 1

    # Syndrome [0,0,0,1,1] is in table (column 4) -> correct qubit 4
    result4 = E.decode!(trues(Hcols), ctd, Bool[0,0,0,1,1])
    @test result4[4] == 1
    @test sum(result4) == 1

    # Missing syndrome returns nothing
    @test isnothing(E.decode!(trues(Hcols), ctd, Bool[0,0,0,0,1]))
end
