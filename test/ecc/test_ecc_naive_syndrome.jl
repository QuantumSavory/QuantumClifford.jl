@testitem "ECC - naive syndrome circuits" tags=[:ecc, :ecc_syndrome_measurement_correctness] begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")
    codes = all_testable_code_instances(; maxn=100)
    include("common/naive_syndrome.jl")
end
