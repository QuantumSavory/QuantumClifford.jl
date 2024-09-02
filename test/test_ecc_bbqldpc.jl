@testitem "Bivariate bicycle (BB) quantum LDPC code" begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC, BBQLDPC

    @testset "Verify number of logical qubits `k` from  Table 3" begin
        # Refer to [bravyi2024high](@cite) for code constructions
        @test code_k(BBQLDPC(9 , 6 , [3 , 1 , 2] , [3 , 1 , 2]))   == 8
        @test code_k(BBQLDPC(15, 3 , [9 , 1 , 2] , [0 , 2 , 7]))   == 8
        @test code_k(BBQLDPC(12, 12, [3 , 2 , 7] , [3 , 1 , 2]))   == 12
        @test code_k(BBQLDPC(12, 6 , [3 , 1 , 2] , [3 , 1 , 2]))   == 12
        @test code_k(BBQLDPC(6 , 6 , [3 , 1 , 2] , [3 , 1 , 2]))   == 12
        @test code_k(BBQLDPC(30, 6 , [9 , 1 , 2] , [3 , 25, 26]))  == 12
        @test code_k(BBQLDPC(21, 18, [3 , 10, 17], [5 , 3 , 19]))  == 16
        @test code_k(BBQLDPC(28, 14, [26, 6 , 8] , [7 , 9 , 20]))  == 24
    end
end
