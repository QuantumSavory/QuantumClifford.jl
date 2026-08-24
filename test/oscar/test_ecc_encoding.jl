@testitem "Oscar ECC encoding circuits" tags=[:ecc, :ecc_encoding, :oscar_required] begin
    include("../ecc/test_ecc_base.jl")
    codes = all_testable_code_instances(; instance_args=oscar_code_instance_args)
    include("../ecc/common/encoding.jl")
end
