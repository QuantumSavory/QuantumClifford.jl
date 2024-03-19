using Test
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC

function is_css_matrix(H)
    for i in 1:size(H, 1)
        for j in 1:size(H, 2)
            e = H[i, j]
                if all(x ∈ [(true, false), (false, true)] for x in H) 
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
    
    # Test Perfect5
    @testset "Steane7" begin
        test_code(Steane7())
    end
end
