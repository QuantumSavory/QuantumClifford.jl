@testitem "ECC code metachecks" tags=[:ecc, :ecc_base] begin
    using QECCore
    using QuantumClifford
    using QuantumClifford
    using QuantumClifford.ECC

    include("test_ecc_base.jl")

    all_codes = all_testable_code_instances()
    function test_all_codes_metachecks(all_codes)
        for c in all_codes
            !isnothing(Mx) && @test iszero(mod.(Mx*parity_matrix_x(c), 2))
            !isnothing(Mz) && @test iszero(mod.(Mz*parity_matrix_z(c), 2))
            isnothing(Mx) && isnothing(Mz) && @info "No metacheck matrices defined for $(typeof(c))"
        end
    end
end
