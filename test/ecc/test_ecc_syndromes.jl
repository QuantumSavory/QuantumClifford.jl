@testitem "ECC Syndromes" tags=[:ecc, :ecc_syndrome_circuit_equivalence] begin
    using QuantumClifford: mul_left!, embed
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")
    codes = all_testable_code_instances()
    include(joinpath(@__DIR__, "common", "syndromes.jl"))
end
