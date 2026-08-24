@testitem "ECC code properties" tags=[:ecc, :ecc_base] begin
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
            _, _, rank = canonicalize!(copy(H), ranks=true)
            # Parity matrices do not encode stabilizer phases.
            canonical_H = canonicalize!(Stabilizer(stab_to_gf2(H)))
            matrix_H = Stabilizer(parity_matrix(code))
            @test nqubits(code) == size(H, 2) == code_n(code)
            @test size(H, 1) == code_s(code)
            @test size(matrix_H) == size(H)
            @test code_s(code) + code_k(code) >= code_n(code) # possibly exist redundant checks
            @test rank <= size(H, 1)
            @test QuantumClifford.stab_looks_good(copy(H), remove_redundant_rows=true)
            @test QuantumClifford.stab_looks_good(matrix_H, remove_redundant_rows=true)
            @test canonical_H == canonicalize!(matrix_H)

            if iscss(code) === true
                Hx = parity_matrix_x(code)
                Hz = parity_matrix_z(code)
                @test size(Hx, 2) == size(Hz, 2) == code_n(code)
                @test size(Hx, 1) + size(Hz, 1) == code_s(code)
                @test iszero(mod.(Hx * LinearAlgebra.transpose(Hz), 2))

                split_H = parity_checks(CSS(Hx, Hz))
                @test QuantumClifford.stab_looks_good(split_H, remove_redundant_rows=true)
                @test canonical_H == canonicalize!(split_H)
            end
        end
    end

    @testset "is degenerate function - test on popular codes" begin
        @test isdegenerate(Shor9()) == true
        @test isdegenerate(Steane7()) == false
        @test isdegenerate(Steane7(), 2) == true
        @test isdegenerate(Bitflip3()) == true
    end
end
