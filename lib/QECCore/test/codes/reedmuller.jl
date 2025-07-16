@testitem "Reed-Muller" begin
    using Test
    using Nemo: echelon_form, matrix, GF
    using QECCore.LinearAlgebra
    using QECCore
    using QECCore: generator

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
                @test rate(ReedMuller(r, m)) == sum(binomial.(m, 0:r)) / 2 ^ m == rate(RecursiveReedMuller(r, m))
                @test code_n(ReedMuller(r, m)) == 2 ^ m == code_n(RecursiveReedMuller(r, m))
                @test code_k(ReedMuller(r, m)) == sum(binomial.(m, 0:r)) == code_k(RecursiveReedMuller(r, m))
                @test distance(ReedMuller(r, m)) == 2 ^ (m - r) == distance(RecursiveReedMuller(r, m))
                H₁ = generator(RecursiveReedMuller(r, m))
                # generator(RecursiveReedMuller(r, m)) is canonically equivalent to the generator(ReedMuller(r, m)) under reduced row echelon form.
                @test echelon_form(matrix(GF(2), Matrix{Int64}(H))) == echelon_form(matrix(GF(2), Matrix{Int64}(H₁)))
                H = parity_matrix(ReedMuller(r, m))
                H₁ = parity_matrix(RecursiveReedMuller(r, m))
                # parity_matrix(RecursiveReedMuller(r, m)) is canonically equivalent to the parity_matrix(ReedMuller(r, m)) under reduced row echelon form.
                @test echelon_form(matrix(GF(2), Matrix{Int64}(H))) == echelon_form(matrix(GF(2), Matrix{Int64}(H₁)))
                # dim(ReedMuller(m - r - 1, m)) = dim(ReedMuller(r, m)^⊥).
                # ReedMuller(m - r - 1, m) = ReedMuller(r, m)^⊥ ∴ parity check matrix (H) of ReedMuller(r, m) is the generator matrix (G) for ReedMuller(m - r - 1, m).
                H₁ = parity_matrix(ReedMuller(m - r - 1, m))
                G₂ = generator(ReedMuller(r, m))
                @test size(H₁) == size(G₂)
                @test echelon_form(matrix(GF(2), Matrix{Int64}(H₁))) == echelon_form(matrix(GF(2), Matrix{Int64}(G₂)))
                G₃ = generator(ReedMuller(m - r - 1, m))
                H₄ = parity_matrix(ReedMuller(r, m))
                @test size(G₃) == size(H₄)
                # dim(RecursiveReedMuller(m - r - 1, m)) = dim(RecursiveReedMuller(r, m)^⊥).
                # RecursiveReedMuller(m - r - 1, m) = RecursiveReedMuller(r, m)^⊥ ∴ parity check matrix (H) of RecursiveReedMuller(r, m) is the generator matrix (G) for RecursiveReedMuller(m - r - 1, m).
                H₁ = parity_matrix(RecursiveReedMuller(m - r - 1, m))
                G₂ = generator(RecursiveReedMuller(r, m))
                @test echelon_form(matrix(GF(2), Matrix{Int64}(H₄))) == echelon_form(matrix(GF(2), Matrix{Int64}(G₃)))
                @test echelon_form(matrix(GF(2), Matrix{Int64}(H₁))) == echelon_form(matrix(GF(2), Matrix{Int64}(G₂)))
                G₃ = generator(RecursiveReedMuller(m - r - 1, m))
                H₄ = parity_matrix(RecursiveReedMuller(r, m))
                @test size(G₃) == size(H₄)
                @test echelon_form(matrix(GF(2), Matrix{Int64}(H₄))) == echelon_form(matrix(GF(2), Matrix{Int64}(G₃)))
            end
        end
    end

    @testset "Test special case 1: generator(RM(0, m))" begin
        for m in 3:10
            H = generator(ReedMuller(0, m))
            expected = ones(Int64, 1, 2^m)
            @test H == expected
        end
    end

    @testset "Test special case 2: generator(RM(m, m))" begin
        for m in 3:10
            H = generator(ReedMuller(m, m))
            expected = check_RM_m_m(m)
            @test echelon_form(matrix(GF(2), Matrix{Int64}(H))) == echelon_form(matrix(GF(2), expected))
        end
    end

    @testset "Test common examples of RM(r,m) codes" begin
        # Examples from [raaphorst2003reed](@cite), [djordjevic2021quantum](@cite), [abbe2020reed](@cite).
        # RM(0,3)
        @test generator(ReedMuller(0,3)) == [1 1 1 1 1 1 1 1]
        # RM(1,3)
        @test generator(ReedMuller(1,3)) == [1 1 1 1 1 1 1 1;
                                             1 1 1 1 0 0 0 0;
                                             1 1 0 0 1 1 0 0;
                                             1 0 1 0 1 0 1 0]
        # RM(2,3)
        @test generator(ReedMuller(2,3)) == [1 1 1 1 1 1 1 1;
                                             1 1 1 1 0 0 0 0;
                                             1 1 0 0 1 1 0 0;
                                             1 0 1 0 1 0 1 0;
                                             1 1 0 0 0 0 0 0;
                                             1 0 1 0 0 0 0 0;
                                             1 0 0 0 1 0 0 0]
        # RM(3,3)
        @test generator(ReedMuller(3,3)) == [1 1 1 1 1 1 1 1;
                                             1 1 1 1 0 0 0 0;
                                             1 1 0 0 1 1 0 0;
                                             1 0 1 0 1 0 1 0;
                                             1 1 0 0 0 0 0 0;
                                             1 0 1 0 0 0 0 0;
                                             1 0 0 0 1 0 0 0;
                                             1 0 0 0 0 0 0 0]
        # RM(2,4)
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
end
