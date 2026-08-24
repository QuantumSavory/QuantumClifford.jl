@testitem "Oscar ECC syndromes" tags=[:ecc, :ecc_syndrome_circuit_equivalence, :oscar_required] begin
    using QuantumClifford: mul_left!, embed
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("../ecc/test_ecc_base.jl")
    codes = testable_code_instances(oscar_code_instance_args)
    include("../ecc/common/syndromes.jl")
end
