using Test
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC

function is_css_matrix(H)
    nrows, ncols = size(H)
    for i in 1:nrows
        has_x = false
        has_z = false
        for j in 1:ncols
            e = H[i, j]
            has_x |= H[i,j][1]
            has_z |= H[i,j][2]
            has_x && has_z && return false
        end
    end
    return true
end

function test_code(code)
    H = parity_checks(code)
    @test iscss(code) == is_css_matrix(H)
    if iscss(code)
        @test is_css_matrix(H) == true "Parity check matrix for CSS code should be in CSS form."
    end
end

known_all_codes = [Shor9(), Steane7(), Gottesman(3), Cleve8(), Perfect5(), Toric(8,8), CSS([0 1 1 0; 1 1 0 0], [1 1 1 1]), Bitflip3()]

@testset "is CSS" begin
    for code in known_all_codes
        test_code(code)
    end
end
