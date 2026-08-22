@testitem "BCH namespace compatibility" tags=[:ecc, :ecc_base] begin
    import QECCore
    import QuantumClifford

    @test QuantumClifford.ECC.BCH === QECCore.BCH
    @test QuantumClifford.ECC.AbstractPolynomialCode === QECCore.AbstractPolynomialCode
    @test QuantumClifford.ECC.BCH(3, 1) isa QECCore.BCH
end
