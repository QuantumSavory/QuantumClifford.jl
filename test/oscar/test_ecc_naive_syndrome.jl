@testitem "Oscar ECC naive syndrome circuits" tags=[:ecc, :ecc_syndrome_measurement_correctness, :oscar_required] begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include(joinpath(@__DIR__, "..", "ecc", "test_ecc_base.jl"))
    codes = testable_code_instances(oscar_code_instance_args; maxn=100)
    include(joinpath(@__DIR__, "..", "ecc", "common", "naive_syndrome.jl"))
end
