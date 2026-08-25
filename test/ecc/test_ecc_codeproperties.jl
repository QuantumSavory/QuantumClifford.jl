@testitem "ECC code properties" tags=[:ecc, :ecc_base] begin
    include("test_ecc_base.jl")
    codes = all_testable_code_instances()
    include("common/codeproperties.jl")
end
