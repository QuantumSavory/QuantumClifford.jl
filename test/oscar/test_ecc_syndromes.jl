@testitem "Oscar ECC syndromes" tags=[:ecc, :ecc_syndrome_circuit_equivalence, :oscar_required] begin
    include("../ecc/test_ecc_base.jl")
    include("test_ecc_oscar_base.jl")
    codes = all_testable_code_instances(; instance_args=oscar_code_instance_args)
    include("../ecc/common/syndromes.jl")
end
