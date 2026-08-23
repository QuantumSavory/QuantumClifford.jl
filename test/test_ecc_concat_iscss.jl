@testitem "Concatenated CSS codes are identified as CSS" tags=[:ecc, :ecc_base] begin
    using QuantumClifford.ECC: Concat, Steane7, iscss

    @test iscss(Concat(Steane7(), Steane7())) === true
end
