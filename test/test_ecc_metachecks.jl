@testitem "ECC code metachecks" tags=[:ecc, :ecc_base] begin
    using QECCore
    using QuantumClifford
    using QuantumClifford
    using QuantumClifford.ECC

    include("test_ecc_base.jl")

    all_codes = all_testable_code_instances()
    
    function test_all_codes_metachecks(all_codes)
        for c in all_codes
            code_mx = (hasproperty(c, :D) && c.D < 4) || # for DDimensional Toric and surface codes, Mx exist for D = 4
                (hasproperty(c, :boundary_maps) && length(c.boundary_maps) != 4) || # for homological product, length 4 means Mx exist
                (hasproperty(c, :t) && c.t < 4) || # for MM code, t  = 4 means Mx exist, otherwise not
                (hasproperty(c, :polynomials) && length(c.polynomials) != 4)

            code_mz = (hasproperty(c, :D) && c.D < 3) || # for DDimensional Toric and surface codes, Mz exist for D = 3
                (hasproperty(c, :boundary_maps) && length(c.boundary_maps) != 3) || # for homological product, length 3 means Mz exist
                (hasproperty(c, :t) && c.t < 3) || # for MM code, 3  = 4 means Mx exist, otherwise not
                (hasproperty(c, :polynomials) && length(c.polynomials) != 4) # for TT code, Mz exist since polynomials are 3

            @testset "Metacheck tests for $(typeof(c))" begin
                if hasmethod(metacheck_matrix_x, Tuple{typeof(c)}) && !code_mx
                    Hx = parity_matrix_x(c)
                    Mx = metacheck_matrix_x(c)
                    @test iszero(mod.(Mx*Hx, 2))
                end
                if hasmethod(metacheck_matrix_z, Tuple{typeof(c)}) && !code_mz
                    Hz = parity_matrix_z(c)
                    Mz = metacheck_matrix_z(c)
                    @test iszero(mod.(Mz*Hz, 2))
                end
            end
        end
    end
end
