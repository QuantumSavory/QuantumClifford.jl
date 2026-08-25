@testitem "encoding circuits - compare to algebraic construction of encoded state" tags=[:ecc, :ecc_encoding] begin
    include("test_ecc_base.jl")
    codes = [
        all_testable_code_instances()...,
        S"Y_",
        S"Z_",
        S"X_",
        [random_stabilizer(5,7) for _ in 1:100]...
    ]
    include("common/encoding.jl")
end
