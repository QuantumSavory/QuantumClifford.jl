@testitem "ECC code properties" tags=[:ecc, :ecc_base] begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")
    codes = all_testable_code_instances()
    include(joinpath(@__DIR__, "common", "codeproperties.jl"))
end
