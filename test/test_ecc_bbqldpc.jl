@testitem "Bivariate bicycle (BB) quantum LDPC code" begin
    using QuantumClifford: stab_to_gf2
    using QuantumClifford.ECC
    using Nemo: nullspace, GF, matrix, hom, free_module, kernel, domain, map, gens
    using QuantumClifford.ECC: AbstractECC, BBQLDPC, parity_checks_x, parity_checks_z

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

    @testset "Verify number of logical qubits `k` from  Table 3" begin
        # Refer to [bravyi2024high](@cite) for code constructions
        @test code_k(BBQLDPC(9 , 6 , [3 , 1 , 2] , [3 , 1 , 2]))   == 8   == _formula_k(BBQLDPC(9 , 6 , [3 , 1 , 2] , [3 , 1 , 2]))
        @test code_k(BBQLDPC(15, 3 , [9 , 1 , 2] , [0 , 2 , 7]))   == 8   == _formula_k(BBQLDPC(15, 3 , [9 , 1 , 2] , [0 , 2 , 7]))
        @test code_k(BBQLDPC(12, 12, [3 , 2 , 7] , [3 , 1 , 2]))   == 12  == _formula_k(BBQLDPC(12, 12, [3 , 2 , 7] , [3 , 1 , 2]))
        @test code_k(BBQLDPC(12, 6 , [3 , 1 , 2] , [3 , 1 , 2]))   == 12  == _formula_k(BBQLDPC(12, 6 , [3 , 1 , 2] , [3 , 1 , 2]))
        @test code_k(BBQLDPC(6 , 6 , [3 , 1 , 2] , [3 , 1 , 2]))   == 12  == _formula_k(BBQLDPC(6 , 6 , [3 , 1 , 2] , [3 , 1 , 2]))
        @test code_k(BBQLDPC(30, 6 , [9 , 1 , 2] , [3 , 25, 26]))  == 12  == _formula_k(BBQLDPC(30, 6 , [9 , 1 , 2] , [3 , 25, 26]))
        @test code_k(BBQLDPC(21, 18, [3 , 10, 17], [5 , 3 , 19]))  == 16  == _formula_k(BBQLDPC(21, 18, [3 , 10, 17], [5 , 3 , 19]))
        @test code_k(BBQLDPC(28, 14, [26, 6 , 8] , [7 , 9 , 20]))  == 24  == _formula_k(BBQLDPC(28, 14, [26, 6 , 8] , [7 , 9 , 20]))
    end
end
