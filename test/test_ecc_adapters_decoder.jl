@testitem "Adapters decoder integration" tags=[:ecc, :ecc_decoding] begin
    # End-to-end sanity check: wrap a merged code in AdapterMergedCode,
    # build a BP-OSD decoder via PyQDecoders, run a CommutationCheckECCSetup
    # Monte Carlo evaluation, assert the logical error rate is finite,
    # in [0, 1], and roughly consistent with a distance-3 code at code-
    # capacity noise. Guards against regressions in the wrapper that would
    # silently break decoder integration.
    using QuantumClifford
    using QuantumClifford.ECC: CSS, Surface, parity_matrix_x, parity_matrix_z,
                               logz_ops, code_n, code_k,
                               PyBeliefPropOSDecoder, evaluate_decoder,
                               CommutationCheckECCSetup
    using QuantumClifford.ECC.Adapters: CodePair, build_adapter, AdapterMergedCode
    using QuantumClifford: stab_to_gf2
    using Random

    import PyQDecoders   # load the extension

    as_css(c) = CSS(Matrix{Bool}(parity_matrix_x(c)),
                    Matrix{Bool}(parity_matrix_z(c)))

    c1 = as_css(Surface(3, 3)); c2 = as_css(Surface(3, 3))
    lz = stab_to_gf2(logz_ops(Surface(3, 3)))
    n0 = code_n(c1)
    z = sort(findall(!iszero, lz[1, n0+1:2n0]))
    adapter = build_adapter(CodePair(c1, c2, z, z))
    wrap = AdapterMergedCode(adapter)

    @test code_n(wrap) == 33
    @test code_k(wrap) == 1

    p = 0.01
    Random.seed!(2026)
    dec = PyBeliefPropOSDecoder(wrap; maxiter = 200, bpmethod = :minsum,
                                 errorrate = p, osdmethod = :exhaustive,
                                 osdorder = 8)
    nshots = 500
    result = evaluate_decoder(dec, CommutationCheckECCSetup(p), nshots)

    @test isa(result, Real)
    @test isfinite(result)
    @test 0.0 ≤ result ≤ 0.05   # d=3 code at p=0.01 → expect p_L ≈ 1e-3..1e-2
end
