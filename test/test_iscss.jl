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
            if e == (true, false)
                has_x = true
            elseif e == (false, true)
                has_z = true
            end
            if has_x && has_z
                return false
            end
        end
    end
    return true
end

function test_code(code)
    is_css = iscss(code)
    H = parity_checks(code)
    is_css_matrix_result = is_css_matrix(H)
    println("Is CSS: ", is_css)
    println("Is parity check matrix in CSS form: ", is_css_matrix_result)
end

@testset "Test code and parity check matrix" begin
    @testset "Shor9" begin
        test_code(Shor9())
    end
    
    # Test Steane7
    @testset "Steane7" begin
        test_code(Steane7())
    end

    # Test Gottesman(3)
    @testset "Gottesman" begin
        test_code(Gottesman(3))
    end
end
