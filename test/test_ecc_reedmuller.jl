using Test
using Nemo: echelon_form, QQ, matrix, GF
using LinearAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedMuller, generator, RecursiveReedMuller

function designed_distance(matrix, m, r)
    distance = 2 ^ (m - r)
    for row in eachrow(matrix)
        count = sum(row)
        if count < distance 
            return false
        end
    end
    return true
end

# Generate binary matrix representing all possible binary strings of length `2ᵐ` encompassing the entire field `F(2ᵐ)`.
function check_RM_m_m(m::Int)
    n = 2^m
    matrix = Matrix{Bool}(undef, n, n)
    for i in 0:(n-1)
        for j in 0:(n-1)
            matrix[i+1, j+1] = ((i & j) == j)
        end
    end
    return matrix
end

@testset "Test RM(r, m) properties" begin
    for m in 3:10
        for r in 3:m-1
            H = generator(ReedMuller(r, m))
            mat = matrix(GF(2), H)
            computed_rank = LinearAlgebra.rank(mat)
            expected_rank = sum(binomial.(m, 0:r))
            @test computed_rank == expected_rank
            @test designed_distance(H, m, r) == true
            @test rate(ReedMuller(r, m)) == sum(binomial.(m, 0:r)) / 2 ^ m
            @test code_n(ReedMuller(r, m)) == 2 ^ m
            @test code_k(ReedMuller(r, m)) == sum(binomial.(m, 0:r))
            @test distance(ReedMuller(r, m)) == 2 ^ (m - r)
            H₁ = generator(RecursiveReedMuller(r, m))
            @test echelon_form(matrix(QQ, Matrix{Int64}(H))) == echelon_form(matrix(QQ, Matrix{Int64}(H₁)))
            H = parity_checks(ReedMuller(r, m))
            H₁ = parity_checks(RecursiveReedMuller(r, m))
            @test echelon_form(matrix(QQ, Matrix{Int64}(H))) == echelon_form(matrix(QQ, Matrix{Int64}(H₁)))
            H₁ = generator(ReedMuller(m - r - 1, m))
            H = parity_checks(ReedMuller(r, m))
            @test size(H) == size(H₁) # dim(RM(m - r - 1, m)) = dim(RM(r, m)^⊥)
        end
    end

    # Test special case 1: generator(RM(0, m))
    for m in 3:10
        H = generator(ReedMuller(0, m))
        expected = ones(Int64, 1, 2^m)
        @test H == expected
    end

    # Test special case 2: generator(RM(m, m))
    for m in 3:10
        H = generator(ReedMuller(m, m))
        expected = check_RM_m_m(m)
        @test echelon_form(matrix(QQ, Matrix{Int64}(H))) == echelon_form(matrix(QQ, expected))
    end
    
    # Testing common examples of RM(r,m) codes [raaphorst2003reed](@cite), [djordjevic2021quantum](@cite), [abbe2020reed](@cite).
    # RM(0,3)  
    @test generator(ReedMuller(0,3)) == [1 1 1 1 1 1 1 1]
    
    #RM(1,3) 
    @test generator(ReedMuller(1,3)) == [1 1 1 1 1 1 1 1;
                                         1 1 1 1 0 0 0 0;
                                         1 1 0 0 1 1 0 0;
                                         1 0 1 0 1 0 1 0]
    #RM(2,3)
    @test generator(ReedMuller(2,3)) == [1 1 1 1 1 1 1 1;
                                         1 1 1 1 0 0 0 0;
                                         1 1 0 0 1 1 0 0;
                                         1 0 1 0 1 0 1 0;
                                         1 1 0 0 0 0 0 0;
                                         1 0 1 0 0 0 0 0;
                                         1 0 0 0 1 0 0 0]
    #RM(3,3)
    @test generator(ReedMuller(3,3)) == [1 1 1 1 1 1 1 1;
                                         1 1 1 1 0 0 0 0;
                                         1 1 0 0 1 1 0 0;
                                         1 0 1 0 1 0 1 0;
                                         1 1 0 0 0 0 0 0;
                                         1 0 1 0 0 0 0 0;
                                         1 0 0 0 1 0 0 0;
                                         1 0 0 0 0 0 0 0]
    #RM(2,4)
    @test generator(ReedMuller(2,4)) == [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
                                         1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
                                         1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0;
                                         1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0;
                                         1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0;
                                         1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
                                         1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
                                         1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0;
                                         1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
                                         1 0 1 0 0 0 0 0 1 0 1 0 0 0 0 0;
                                         1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0]
end
