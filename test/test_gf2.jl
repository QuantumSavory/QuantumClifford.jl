@testitem "GF(2) representations" begin
    using Random
    using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    @testset "Equivalence of GF(2) Gaussian elimination and Stabilizer canonicalization" begin
        for n in test_sizes
            for rep in 1:5
                s = random_stabilizer(n)[randperm(n)[1:rand(n÷2+1:n)]]
                cs = canonicalize!(copy(s));
                H = stab_to_gf2(cs);
                cH = gf2_gausselim!(stab_to_gf2(s));
                @test H==cH
            end
        end
    end
    @testset "GF(2) H and G matrices" begin
        for n in test_sizes
            for rep in 1:5
                H = random_invertible_gf2(n)[randperm(n)[1:rand(n÷2+1:n)],:]
                H = gf2_gausselim!(H)
                G = gf2_H_to_G(H)
                @test sum(G*H' .%2)==0;
            end
        end
    end
    @testset "GF(2) row echelon form with transformation matrix, pivots etc." begin
        for n in test_sizes
            for rep in 1:100
                gf2_matrices = [rand(Bool, size, size) for size in test_sizes]
                for (i, matrix) in enumerate(gf2_matrices)
                     echelon_form, _, transformation, _ = gf2_row_echelon_with_pivots!(Matrix{Int}(matrix), full=true) # in-place
                     # Check the correctness of the transformation matrix
                     @test (transformation*matrix) .%2 == echelon_form
                     # Check the correctness of Gaussian elimination
                     @test echelon_form == gf2_gausselim!(matrix)
                end
            end
        end
    end
    function is_in_nullspace(A, x)
        # Ensure x is the correct orientation
        if size(x, 1) != size(A, 2)
            x = transpose(x)
        end
        # Perform modulo 2 arithmetic: A * x must be zero mod 2
        if size(x, 2) == 1  # x is a single column vector
            result = A * x
            return all(result .% 2 .== 0)  # Check if A * x = 0 mod 2
        else  # x is a matrix, check each column vector
            for i in 1:size(x, 2)
                result = A * x[:, i] # Multiply A with the i-th column of x
                if !all(result .% 2 .== 0) # Check if A * column = 0 mod 2
                    return false
                end
            end
            return true # All columns are in the null space mod 2
        end
    end
    @testset "GF(2) nullspace of the binary matrix" begin
        for n in test_sizes
            for rep in 1:100
                gf2_matrices = [rand(Bool, size, size) for size in test_sizes]
                for (i, matrix) in enumerate(gf2_matrices)
                    imat = Matrix{Int}(matrix)
                    ns = gf2_nullspace(imat)
                    @test is_in_nullspace(imat, ns)
                end
            end
        end
    end
    @testset "Consistency check with ldpc" begin
        # sanity checks for comparison to https://github.com/quantumgizmos/ldpc
        # results compared with 'from ldpc.mod2 import nullspace, row_basis, row_echelon'
        # Consistency check 1
        H = [1 1 1; 1 1 1; 0 1 0]
        echelon_form, rank, transformation, pivots = gf2_row_echelon_with_pivots!(copy(H)) # in-place
        @test echelon_form == [1 1 1; 0 1 0; 0 0 0]
        @test rank == 2
        @test transformation == [1 0 0; 0 0 1; 1 1 0]
        @test pivots == [1, 2] # in python, it's [0, 1] due to zero-indexing
        @test mod.((transformation*copy(H)), 2) == echelon_form
        @test gf2_nullspace(copy(H)) == [1 0 1]
        @test gf2_rowspace_basis(copy(H)) == [1 1 1; 0 1 0]
        # Consistency check 2
        H = [0 0 0 1 1 1 1;
             0 1 1 0 0 1 1;
             1 0 1 0 1 0 1]
        echelon_form, rank, transformation, pivots = gf2_row_echelon_with_pivots!(copy(H)) # in-place
        @test echelon_form == [1 0 1 0 1 0 1;
                               0 1 1 0 0 1 1;
                               0 0 0 1 1 1 1]
        @test rank == 3
        @test transformation == [0 0 1;
                                 0 1 0;
                                 1 0 0]
        @test pivots == [1, 2, 4] # in python, it's [0, 1, 3] due to zero-indexing
        @test mod.((transformation*copy(H)), 2) == echelon_form
        @test gf2_nullspace(copy(H)) == [1 1 1 0 0 0 0;
                                         0 1 1 1 1 0 0;
                                         0 1 0 1 0 1 0;
                                         0 0 1 1 0 0 1]
        @test gf2_rowspace_basis(copy(H)) == [0 0 0 1 1 1 1;
                                              0 1 1 0 0 1 1;
                                              1 0 1 0 1 0 1]
        # Consistency check 3
        H = [1 1 0; 0 1 1; 1 0 1]
        echelon_form, rank, transformation, pivots = gf2_row_echelon_with_pivots!(copy(H)) # in-place
        @test echelon_form == [1 1 0;
                               0 1 1;
                               0 0 0]
        @test rank == 2
        @test transformation == [1 0 0;
                                 0 1 0;
                                 1 1 1]
        @test pivots == [1,2 ] # in python, it's [0, 1] due to zero-indexing
        @test mod.((transformation*copy(H)), 2) == echelon_form
        @test gf2_nullspace(copy(H)) == [1 1 1]
        @test gf2_rowspace_basis(copy(H)) == [1 1 0;
                                              0 1 1]
        # Consistency check 4
        H = [1 1 0; 0 1 0; 0 0 1]
        echelon_form, rank, transformation, pivots = gf2_row_echelon_with_pivots!(copy(H)) # in-place
        @test echelon_form == [1 1 0;
                               0 1 0;
                               0 0 1]
        @test rank == 3
        @test transformation == [1 0 0;
                                 0 1 0;
                                 0 0 1]
        @test pivots == [1, 2, 3]  # in python, it's [0, 1, 2] due to zero-indexing
        @test mod.((transformation*copy(H)), 2) == echelon_form
        @test gf2_nullspace(copy(H)) == [0 0 0]
        @test gf2_rowspace_basis(copy(H)) == [1 1 0;
                                              0 1 0;
                                              0 0 1]
        # Consistency check 5
        H = [1 1 0; 0 1 0; 0 0 1; 0 1 1]
        echelon_form, rank, transformation, pivots = gf2_row_echelon_with_pivots!(copy(H)) # in-place
        @test echelon_form == [1 1 0;
                               0 1 0;
                               0 0 1;
                               0 0 0]
        @test rank == 3
        @test transformation == [1 0 0 0;
                                 0 1 0 0;
                                 0 0 1 0;
                                 0 1 1 1]
        @test pivots == [1, 2, 3] # in python, it's [0, 1, 2] due to zero-indexing
        @test mod.((transformation*copy(H)), 2) == echelon_form
        @test gf2_nullspace(copy(H)) == [0 0 0]
        @test gf2_rowspace_basis(copy(H)) == [1 1 0;
                                              0 1 0;
                                              0 0 1]
        # Consistency check 6
        H = [0 0 0 1 1 1 1;
             0 1 1 0 0 1 1;
             1 0 1 0 1 0 1]
        echelon_form, rank, transformation, pivots = gf2_row_echelon_with_pivots!(copy(H)) # in-place
        @test echelon_form == [1 0 1 0 1 0 1;
                               0 1 1 0 0 1 1;
                               0 0 0 1 1 1 1]
        @test rank == 3
        @test transformation == [0 0 1;
                                 0 1 0;
                                 1 0 0]
        @test pivots == [1, 2, 4]  # in python, it's [0, 1, 3] due to zero-indexing
        @test mod.((transformation*copy(H)), 2) == echelon_form
        @test gf2_nullspace(copy(H)) == [1 1 1 0 0 0 0;
                                         0 1 1 1 1 0 0;
                                         0 1 0 1 0 1 0;
                                         0 0 1 1 0 0 1]
        @test gf2_rowspace_basis(copy(H)) == [0 0 0 1 1 1 1;
                                              0 1 1 0 0 1 1;
                                              1 0 1 0 1 0 1]
    end
end
