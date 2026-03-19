@testset "Flag ancillary Pauli measurement" begin
    using QuantumClifford

    # simple Pauli operator (example)
    p = random_pauli(5)

    circ, anc, bits = flag_ancillary_paulimeasurement(p)

    @test length(circ) > 0
    @test anc == 2
    @test length(bits) == 2
end
