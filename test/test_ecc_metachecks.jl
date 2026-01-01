@testitem "ECC code metachecks" tags=[:ecc, :ecc_base] begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")

    all_codes = all_testable_code_instances()
    
    function test_all_codes_metachecks(all_codes)
        for c in all_codes
            @testset "Metacheck tests for $(typeof(c))" begin
                try
                    Hx = parity_matrix_x(c)
                    Mx = metacheck_matrix_x(c)
                    @test iszero(mod.(Mx*Hx, 2))
                catch e
                    # X-metacheck not defined for this code instance, e.g. 2D Toric code
                end
                try
                    Hz = parity_matrix_z(c)
                    Mz = metacheck_matrix_z(c)
                    @test iszero(mod.(Mz*Hz, 2))
                catch e
                    # Z-metacheck not defined for this code instance, e.g. 2D Toric code
                end
            end
        end
    end
end
