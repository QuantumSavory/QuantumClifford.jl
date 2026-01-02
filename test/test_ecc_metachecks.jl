@testitem "ECC code metachecks" tags=[:ecc, :ecc_base] begin
    using QECCore
    using QECCore: hasmetachecks
    using QuantumClifford
    using QuantumClifford
    using QuantumClifford.ECC

    include("test_ecc_base.jl")

    all_codes = all_testable_code_instances()

    function test_all_codes_metachecks(all_codes)
        for c in all_codes
            has_mx, has_mz = hasmetachecks(c)
            if has_mz && hasmethod(metacheck_matrix_z, Tuple{typeof(c)})
                Mz = metacheck_matrix_z(c)
                !isempty(Mz) && @test iszero(mod.(Mz*parity_matrix_z(c), 2))
            end
            if has_mx && hasmethod(metacheck_matrix_x, Tuple{typeof(c)})
                Mx = metacheck_matrix_x(c)
                !isempty(Mx) && @test iszero(mod.(Mx*parity_matrix_x(c), 2))
            end
        end
    end
end
