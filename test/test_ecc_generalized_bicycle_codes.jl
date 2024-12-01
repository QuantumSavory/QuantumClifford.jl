@testitem "ECC GB" begin
    using Hecke
    using JuMP
    using GLPK
    using QuantumClifford.ECC: generalized_bicycle_codes, code_k, code_n, distance

    # codes taken from Table 1 of [lin2024quantum](@cite)
    # Abelian 2BGA codes can be viewed as GB codes.

    codes = [generalized_bicycle_codes([0 ,  1,  3,  7], [0 ,  1, 12, 19], 27), # [[54, 6,  9]]
             generalized_bicycle_codes([0 , 10,  6, 13], [0 , 25, 16, 12], 30), # [[60, 6, 10]]
             generalized_bicycle_codes([0 , 15, 16, 18], [0 ,  1, 24, 27], 35), # [[70, 8, 10]]
             generalized_bicycle_codes([0 ,  9, 28, 31], [0 ,  1, 21, 34], 36), # [[72, 8, 10]]
             generalized_bicycle_codes([0 ,  9, 28, 13], [0 ,  1,  3, 22], 36)] # [[72, 8,  9]]

    n_values = [54, 60, 70, 72, 72]
    k_values = [6 , 6 , 8 , 8 , 10]
    d_values = [9 , 10, 10, 10,  9]

    for i in 1:length(n_values)
        @test code_n(codes[i]) == n_values[i]
    end

    for i in 1:length(n_values)
        @test code_k(codes[i]) == k_values[i]
    end

    for i in 1:length(d_values)
        j = rand(1:code_k(codes[i]))
        @test distance(codes[i], logical_qubit=j) == d_values[i]
    end
end
