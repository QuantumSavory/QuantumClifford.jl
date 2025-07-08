@testitem "ECC circulant_bivariate_bicycle" begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Test
        using QuantumClifford: stab_to_gf2
        using QuantumClifford.ECC
        using Nemo: nullspace, GF, matrix
        using Oscar: hom, free_module, kernel, domain, map, gens
        using QuantumClifford.ECC: AbstractECC, circulant_bivariate_bicycle, parity_checks_x, parity_checks_z
        # According to Lemma 1 from [bravyi2024high](@cite), k = 2·dim(ker(A)∩ker(B)).
        function _formula_k(stab)
            Hx = parity_checks_x(stab)
            n = size(Hx,2)÷2
            A = matrix(GF(2), Hx[:,1:n])
            B = matrix(GF(2), Hx[:,n+1:end])
            k = GF(2)
            hA = hom(free_module(k, size(A, 1)), free_module(k, size(A, 2)), A)
            hB = hom(free_module(k, size(B, 1)), free_module(k, size(B, 2)), B)
            ans = kernel(hA)[1] ∩ kernel(hB)[1]
            k = 2*size(map(domain(hA), gens(ans[1])), 1)
            return k
        end

        @testset "Verify number of logical qubits `k` from Table 3: bravyi2024high" begin
            # Test cases from [bravyi2024high](@cite)
            code = circulant_bivariate_bicycle(9, 6, [3, 1, 2], [3, 1, 2])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(15, 3, [9, 1, 2], [0, 2, 7])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(12, 12, [3, 2, 7], [3, 1, 2])
            @test code_k(code) == 12 == _formula_k(code)
            code = circulant_bivariate_bicycle(12, 6, [3, 1, 2], [3, 1, 2])
            @test code_k(code) == 12 == _formula_k(code)
            code = circulant_bivariate_bicycle(6, 6, [3, 1, 2], [3, 1, 2])
            @test code_k(code) == 12 == _formula_k(code)
            code = circulant_bivariate_bicycle(30, 6, [9, 1, 2], [3, 25, 26])
            @test code_k(code) == 12 == _formula_k(code)
            code = circulant_bivariate_bicycle(21, 18, [3, 10, 17], [5, 3, 19])
            @test code_k(code) == 16 == _formula_k(code)
            code = circulant_bivariate_bicycle(28, 14, [26, 6, 8], [7, 9, 20])
            @test code_k(code) == 24 == _formula_k(code)
        end

        @testset "Verify number of logical qubits `k` from Table 1: berthusen2024toward" begin
            # Test cases from [berthusen2024toward](@cite)
            code = circulant_bivariate_bicycle(12, 3, [9, 1, 2], [0, 1, 11])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(9, 5, [8, 4, 1], [5, 8, 7])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(12, 5, [10, 4, 1], [0, 1, 2])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(15, 5, [5, 2, 3], [2, 7, 6])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(14, 7, [6, 5, 6], [0, 4, 13])
            @test code_k(code) == 12 == _formula_k(code)
        end

        @testset "Verify number of logical qubits `k` from Table 1: wang2024coprime" begin
            # Test cases from [wang2024coprime](@cite)
            code = circulant_bivariate_bicycle(3, 9, [0, 2, 4], [3, 1, 2])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(7, 7, [3, 5, 6], [2, 3, 5])
            @test code_k(code) == 6 == _formula_k(code)
            code = circulant_bivariate_bicycle(3, 21, [0, 2, 10], [3, 1, 2])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(5, 15, [0, 6, 8], [5, 1, 4])
            @test code_k(code) == 16 == _formula_k(code)
            code = circulant_bivariate_bicycle(3, 27, [0, 10, 14], [12, 1, 2])
            @test code_k(code) == 8 == _formula_k(code)
            code = circulant_bivariate_bicycle(6, 15, [3, 1, 2], [6, 4, 5])
            @test code_k(code) == 8 == _formula_k(code)
        end
    end
end
