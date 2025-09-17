@testitem "ECC code properties" tags=[:ecc, :ecc_code_properties] begin
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
        for code in all_testable_code_instances()
            H = parity_checks(code)
            @test iscss(code) in (is_css_matrix(H), nothing)
        end
    end

    @testset "code tableau consistency" begin
        for code in all_testable_code_instances()
            H = parity_checks(code)
            @test nqubits(code) == size(H, 2) == code_n(code)
            @test size(H, 1) == code_s(code)
            @test code_s(code) + code_k(code) >= code_n(code) # possibly exist redundant checks
            _, _, rank = canonicalize!(copy(H), ranks=true)
            @test rank <= size(H, 1)
            @test QuantumClifford.stab_looks_good(copy(H), remove_redundant_rows=true)
        end
    end
end
