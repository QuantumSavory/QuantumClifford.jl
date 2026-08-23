@testitem "Python decoders support concurrent reuse" tags=[:ecc, :ecc_decoding] begin
    using QuantumClifford
    using QuantumClifford.ECC
    import PyQDecoders

    code = Steane7()
    decoders = Vector{Any}(undef, 4)
    @sync for i in eachindex(decoders)
        Threads.@spawn decoders[i] = PyMatchingDecoder(code)
    end

    decoder = first(decoders)
    syndrome = falses(size(parity_checks(code), 1))
    expected = decode(decoder, syndrome)
    results = Vector{Any}(undef, 32)
    @sync for i in eachindex(results)
        Threads.@spawn results[i] = decode(decoder, syndrome)
    end
    @test all(==(expected), results)
end

@testitem "Tesseract decoders support concurrent reuse" tags=[:ecc, :ecc_decoding, :tesseract_required] begin
    using QuantumClifford
    using QuantumClifford.ECC
    import PyTesseractDecoder

    code = Steane7()
    decoders = Vector{Any}(undef, 4)
    @sync for i in eachindex(decoders)
        Threads.@spawn decoders[i] = TesseractDecoder(code)
    end

    decoder = first(decoders)
    syndrome = falses(size(parity_checks(code), 1))
    expected = decode(decoder, syndrome)
    results = Vector{Any}(undef, 32)
    @sync for i in eachindex(results)
        Threads.@spawn results[i] = decode(decoder, syndrome)
    end
    @test all(==(expected), results)
end
