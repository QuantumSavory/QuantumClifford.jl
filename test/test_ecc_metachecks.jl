@testitem "ECC code metachecks" tags=[:ecc, :ecc_base] begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")

    all_codes = all_testable_code_instances()
    
    function test_all_codes_metachecks(all_codes)
        for c in all_codes
            @testset "Metacheck tests for $(typeof(c))" begin
                if hasmethod(metacheck_matrix_x, Tuple{typeof(c)})
                    if (hasproperty(c, :D) && c.D ≥ 4) || (hasproperty(c, :boundary_maps) && length(c.boundary_maps) == 4) || (!hasproperty(c, :D) && !hasproperty(c, :boundary_maps))
                        Hx = parity_matrix_x(c)
                        Mx = metacheck_matrix_x(c)
                        @test iszero(mod.(Mx * Hx, 2))
                    end
                end
                if hasmethod(metacheck_matrix_z, Tuple{typeof(c)})
                    if (hasproperty(c, :D) && c.D ≥ 3) || (hasproperty(c, :boundary_maps) && length(c.boundary_map) == 3) || (!hasproperty(c, :D) && !hasproperty(c, :seed_codes))
                        Hz = parity_matrix_z(c)
                        Mz = metacheck_matrix_z(c)
                        @test iszero(mod.(Mz * Hz, 2))
                    end
                end
            end
        end
    end
end
