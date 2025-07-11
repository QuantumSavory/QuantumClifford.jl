@testitem "ECC Hamming" begin

    using QECCore.LinearAlgebra
    using QECCore
    using Nemo: matrix, GF, echelon_form
    using Test

    function minimum_distance(H)
        n = size(H, 2)
        min_dist = n + 1
        for x_bits in 1:(2^n - 1)
            x = reverse(digits(x_bits, base=2, pad=n))
            xᵀ = reshape(x, :, 1)
            if all(mod.(H * xᵀ, 2) .== 0)
                weight = sum(x)
                min_dist = min(min_dist, weight)
            end
        end
        return min_dist
    end

    @testset "Testing Hamming codes properties" begin
        for r in 3:20
            n = 2 ^ r - 1
            k = 2 ^ r - 1 - r
            H = parity_matrix(Hamming(r))
            H = Matrix{Bool}(H)
            mat = matrix(GF(2), parity_matrix(Hamming(r)))
            computed_rank = rank(mat)
            @test computed_rank == n - k
        end
        # minimum distance test is expensive (NP-hard) so we test for r = [3,4]
        r_vals = [3,4]
        for r in r_vals
            H = parity_matrix(Hamming(r))
            H = Matrix{Bool}(H)
            @test minimum_distance(H) == 3
        end
        # Example taken from [huffman2010fundamentals](@cite).
        @test Matrix{Bool}(parity_matrix(Hamming(3))) == [0  0  0  1  1  1  1;
                                                          0  1  1  0  0  1  1;
                                                          1  0  1  0  1  0  1]
    end
end
