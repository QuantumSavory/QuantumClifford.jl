@testitem "ECC Golay" begin

    using LinearAlgebra
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC, Golay, generator
    using Nemo: matrix, GF, echelon_form

    # Theorem: Let `C` be a binary linear code. If `C` is self-orthogonal and
    # has a generator matrix `G` where each row has weight divisible by four,
    # then every codeword of `C` has weight divisible by four. `H₂₄` is self-dual
    # because its generator matrix has all rows with weight divisible by four.
    # Thus, all codewords of `H₂₄` must have weights divisible by four. Refer to
    # pg. 30 to 33 of Ch1 of Fundamentals of Error Correcting Codes by Huffman,
    # Cary and Pless, Vera.
    function code_weight_property(matrix)
        for row in eachrow(matrix)
            count = sum(row)
            if count % 4 == 0
                return true
            end
        end
        return false
    end

    # Test the equivalence of punctured and extended code by verifying that puncturing
    # the binary parity check matrix H₂₄ in any coordinate and then extending it by adding
    # an overall parity check in the same position yields the original matrix H₂₄.
    # Steps:
    # 1) Puncture the Code: Remove the i-th column from H₂₄ to create a punctured matrix
    # H₂₃. Note: H₂₃ = H[:, [1:i-1; i+1:end]]
    # 2) Extend the Code: Add a column in the same position to ensure each row has even
    # parity. Note: H'₂₄ = [H₂₃[:, 1:i-1] c H₂₃[:, i:end]]. Here, c is a column vector
    # added to ensure each row in H'₂₄ has even parity.
    # 3) Equivalence Check: Verify that H'₂₄ = H₂₄.
    function puncture_code(mat, i)
        return mat[:, [1:i-1; i+1:end]]
    end

    function extend_code(mat, i)
        k, _ = size(mat)
        extended_mat = hcat(mat[:, 1:i-1], zeros(Bool, k, 1), mat[:, i:end])
        # Calculate the parity for each row
        for row in 1:k
            row_parity = sum(mat[row, :]) % 2 == 1
            extended_mat[row, i] = row_parity
        end
        return extended_mat
    end

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

    @testset "Testing binary Golay codes properties" begin
        test_cases = [(24, 12), (23, 12)]
        for (n, k) in test_cases
            H = parity_checks(Golay(n))
            mat = matrix(GF(2), parity_checks(Golay(n)))
            computed_rank = rank(mat)
            @test computed_rank == n - k
        end

        # [24, 12, 8] binary Golay code is a self-dual code [huffman2010fundamentals](@cite).
        H = parity_checks(Golay(24))
        @test code_weight_property(H) == true
        @test H[:, (12 + 1):end]  == H[:, (12 + 1):end]'
        # Example taken from [huffman2010fundamentals](@cite).
        @test parity_checks(Golay(24)) == [0  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0;
                                           1  1  1  0  1  1  1  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  0  0;
                                           1  1  0  1  1  1  0  0  0  1  0  1  0  0  1  0  0  0  0  0  0  0  0  0;
                                           1  0  1  1  1  0  0  0  1  0  1  1  0  0  0  1  0  0  0  0  0  0  0  0;
                                           1  1  1  1  0  0  0  1  0  1  1  0  0  0  0  0  1  0  0  0  0  0  0  0;
                                           1  1  1  0  0  0  1  0  1  1  0  1  0  0  0  0  0  1  0  0  0  0  0  0;
                                           1  1  0  0  0  1  0  1  1  0  1  1  0  0  0  0  0  0  1  0  0  0  0  0;
                                           1  0  0  0  1  0  1  1  0  1  1  1  0  0  0  0  0  0  0  1  0  0  0  0;
                                           1  0  0  1  0  1  1  0  1  1  1  0  0  0  0  0  0  0  0  0  1  0  0  0;
                                           1  0  1  0  1  1  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0;
                                           1  1  0  1  1  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0;
                                           1  0  1  1  0  1  1  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1]

        # minimum distance test
        # [24, 12, 8]
        H = parity_checks(Golay(24))
        @test minimum_distance(H) == 8
        # [23, 12, 7]
        H = parity_checks(Golay(23))
        @test minimum_distance(H) == 7

        # cross-verifying the canonical equivalence of bordered reverse circulant matrix (A)
        # from [huffman2010fundamentals](@cite) with matrix A taken from [bhatia2018mceliece](@cite).
        A = [1 1 0 1 1 1 0 0 0 1 0 1;
             1 0 1 1 1 0 0 0 1 0 1 1;
             0 1 1 1 0 0 0 1 0 1 1 1;
             1 1 1 0 0 0 1 0 1 1 0 1;
             1 1 0 0 0 1 0 1 1 0 1 1;
             1 0 0 0 1 0 1 1 0 1 1 1;
             0 0 0 1 0 1 1 0 1 1 1 1;
             0 0 1 0 1 1 0 1 1 1 0 1;
             0 1 0 1 1 0 1 1 1 0 0 1;
             1 0 1 1 0 1 1 1 0 0 0 1;
             0 1 1 0 1 1 1 0 0 0 1 1;
             1 1 1 1 1 1 1 1 1 1 1 0]

        H = parity_checks(Golay(24))
        @test echelon_form(matrix(GF(2), A)) == echelon_form(matrix(GF(2), H[1:12, 1:12]))
        # test self-duality for extended Golay code, G == H
        @test echelon_form(matrix(GF(2), generator(Golay(24)))) == echelon_form(matrix(GF(2), parity_checks(Golay(24))))

        # All punctured and extended matrices are equivalent to the H₂₄. Test each column for puncturing and extending.
        for i in 1:24
            punctured_mat = puncture_code(H, i)
            extended_mat = extend_code(punctured_mat, i)
            @test extended_mat == H
        end
    end
end
