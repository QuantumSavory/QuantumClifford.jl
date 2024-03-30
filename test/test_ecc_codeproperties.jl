using Test
using QuantumClifford
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
        @test iscss(code) == is_css_matrix(H)
    end
end

@testset "code tableau consistency" begin
    for code in all_testablable_code_instances()
        H = parity_checks(code)
        @test nqubits(code) == size(H, 2) == code_n(code)
        @test size(H, 1) == code_s(code)
        @test code_s(code) + code_k(code) == code_n(code)
        @test size(H, 1) < size(H, 2)
        _, _, rank = canonicalize!(copy(H), ranks=true)
        @test rank == size(H, 1) # TODO maybe weaken this if we want to permit codes with redundancies
        @test QuantumClifford.stab_looks_good(copy(H))
    end
end
