@testitem "Oscar ECC naive syndrome circuits" tags=[:ecc, :ecc_syndrome_measurement_correctness, :oscar_required] begin
    include("../ecc/test_ecc_base.jl")
    codes = all_testable_code_instances(; maxn=100, instance_args=oscar_code_instance_args)
    include("../ecc/common/naive_syndrome.jl")
end
