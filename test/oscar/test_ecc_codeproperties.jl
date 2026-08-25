@testitem "Oscar ECC code properties" tags=[:ecc, :ecc_base, :oscar_required] begin
    include("../ecc/test_ecc_base.jl")
    codes = all_testable_code_instances(; instance_args=oscar_code_instance_args)
    include("../ecc/common/codeproperties.jl")
end
