@testitem "Oscar ECC code properties" tags=[:ecc, :ecc_base, :oscar_required] begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include(joinpath(@__DIR__, "..", "ecc", "test_ecc_base.jl"))
    codes = testable_code_instances(oscar_code_instance_args)
    include(joinpath(@__DIR__, "..", "ecc", "common", "codeproperties.jl"))
end
