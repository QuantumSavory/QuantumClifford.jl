@testitem "ECC belief-propagation decoder" tags=[:ecc, :ecc_decoding, :too_slow_for_ci] begin
    using QuantumClifford.ECC

    import LDPCDecoders
    import PyQDecoders

    include("../ecc/test_ecc_base.jl")

    codes = vcat(LP04, LP118, test_gb_codes, other_lifted_product_codes)
    noise = 0.001
    setups = [
        CommutationCheckECCSetup(noise),
        NaiveSyndromeECCSetup(noise, 0),
        ShorSyndromeECCSetup(noise, 0),
    ]

    for code in codes, setup in setups
        decoder = PyBeliefPropOSDecoder(code; maxiter=2)
        errors = evaluate_decoder(decoder, setup, 10_000)
        @test max(errors...) <= noise
    end
end
