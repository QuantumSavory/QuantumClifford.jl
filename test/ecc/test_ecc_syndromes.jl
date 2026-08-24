@testitem "ECC Syndromes" tags=[:ecc, :ecc_syndrome_circuit_equivalence] begin
    include("test_ecc_base.jl")
    codes = all_testable_code_instances()
    include("common/syndromes.jl")
end
