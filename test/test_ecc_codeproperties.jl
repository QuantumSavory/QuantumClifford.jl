@testitem "ECC code properties" begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    include("test_ecc_base.jl")

    function is_css_matrix(H)
        nrows, ncols = size(H)
        for i in 1:nrows
            has_x = false
            has_z = false
            for j in 1:ncols
                has_x |= H[i,j][1]
                has_z |= H[i,j][2]
                has_x && has_z && return false
            end
        end
        return true
    end

    @testset "is CSS" begin
        for code in all_testablable_code_instances()
            H = parity_checks(code)
            @test iscss(code) in (is_css_matrix(H), nothing)
        end
    end

    @testset "code tableau consistency" begin
        for code in all_testablable_code_instances()
            H = parity_checks(code)
            @test nqubits(code) == size(H, 2) == code_n(code)
            @test size(H, 1) == code_s(code)
            @test code_s(code) + code_k(code) >= code_n(code) # possibly exist redundant checks
            @test size(H, 1) <= size(H, 2) # equal in some cases with redundancy
            @test QuantumClifford.stab_looks_good(copy(H))
        end
    end
end
