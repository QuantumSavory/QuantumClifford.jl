@testitem "Oscar ECC encoding circuits" tags=[:ecc, :ecc_encoding, :oscar_required] begin
    using QuantumClifford
    using QuantumClifford.ECC

    include(joinpath(@__DIR__, "..", "ecc", "test_ecc_base.jl"))
    codes = testable_code_instances(oscar_code_instance_args)
    include(joinpath(@__DIR__, "..", "ecc", "common", "encoding.jl"))
end
