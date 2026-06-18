@testitem "Adapters decoder integration" tags=[:ecc, :ecc_decoding] begin
    # End-to-end sanity check: the merged Adapter plugs directly into the
    # BP-OSD decoder via PyQDecoders, run a CommutationCheckECCSetup Monte
    # Carlo evaluation, assert the logical error rate is finite, in [0, 1],
    # and roughly consistent with a distance-3 code at code-capacity noise.
    using QuantumClifford
    using QuantumClifford.ECC: Surface, logz_ops, code_n, code_k,
                               build_adapter,
                               PyBeliefPropOSDecoder, evaluate_decoder,
                               CommutationCheckECCSetup
    using Random

    import PyQDecoders   # load the extension

    z = logz_ops(Surface(3, 3))[1]
    adapter = build_adapter(Surface(3, 3), Surface(3, 3), z, z)

    @test code_n(adapter) == 33
    @test code_k(adapter) == 1

    p = 0.01
    Random.seed!(2026)
    dec = PyBeliefPropOSDecoder(adapter; maxiter = 200, bpmethod = :minsum,
                                 errorrate = p, osdmethod = :exhaustive,
                                 osdorder = 8)
    nshots = 500
    result = evaluate_decoder(dec, CommutationCheckECCSetup(p), nshots)

    @test isa(result, Real)
    @test isfinite(result)
    @test 0.0 ≤ result ≤ 0.05   # d=3 code at p=0.01 → expect p_L ≈ 1e-3..1e-2
end
