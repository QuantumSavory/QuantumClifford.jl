@testitem "ECC 2BGA" begin
    using Hecke
    using QuantumClifford.ECC: generalized_bicycle_codes, code_k, code_n

    # codes taken from Table 1 of [lin2024quantum](@cite)
    @test code_n(generalized_bicycle_codes([0 , 15, 16, 18], [0 ,  1, 24, 27], 35)) == 70
    @test code_k(generalized_bicycle_codes([0 , 15, 16, 18], [0 ,  1, 24, 27], 35)) == 8
    @test code_n(generalized_bicycle_codes([0 ,  1,  3,  7], [0 ,  1, 12, 19], 27)) == 54
    @test code_k(generalized_bicycle_codes([0 ,  1,  3,  7], [0 ,  1, 12, 19], 27)) == 6
    @test code_n(generalized_bicycle_codes([0 , 10,  6, 13], [0 , 25, 16, 12], 30)) == 60
    @test code_k(generalized_bicycle_codes([0 , 10,  6, 13], [0 , 25, 16, 12], 30)) == 6
    @test code_n(generalized_bicycle_codes([0 ,  9, 28, 31], [0 ,  1, 21, 34], 36)) == 72
    @test code_k(generalized_bicycle_codes([0 ,  9, 28, 31], [0 ,  1, 21, 34], 36)) == 8
    @test code_n(generalized_bicycle_codes([0 ,  9, 28, 13], [0 ,  1, 21, 34], 36)) == 72
    @test code_k(generalized_bicycle_codes([0 ,  9, 28, 13], [0 ,  1,  3, 22], 36)) == 10
end
