using Test
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC

function is_css_matrix(H)
    for row in eachrow(H)
        if all(x ∈ [(true, false), (false, true)] for x in row)
            return false
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
        test_code(Shor9)
    end
    
    # Test Perfect5
    @testset "Perfect5" begin
        test_code(Perfect5)
    end
end
