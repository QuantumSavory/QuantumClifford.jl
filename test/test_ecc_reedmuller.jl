using Test
using Nemo
using Combinatorics
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, ReedMuller

function binomial_coeff_sum(r, m)
    total = 0
    for i in 0:r
        total += length(collect(combinations(1:m, i)))
    end
    return total
end

@testset "Test RM(r, m) Matrix Rank" begin
    for m in 2:5
        for r in 0:m - 1
            H = parity_checks(ReedMuller(r, m))
            mat = Nemo.matrix(Nemo.GF(2), H)
            computed_rank = LinearAlgebra.rank(mat)
            expected_rank = binomial_coeff_sum(r, m)
            @test computed_rank == expected_rank
        end
    end
end

@testset "Testing common examples of RM(r,m) codes [raaphorst2003reed](@cite), [djordjevic2021quantum](@cite), [abbe2020reed](@cite)" begin
    
    #RM(0,3)  
    @test parity_checks(ReedMuller(0,3)) == [1 1 1 1 1 1 1 1]
    
    #RM(1,3) 
    @test parity_checks(ReedMuller(1,3)) == [1 1 1 1 1 1 1 1;
                                             1 1 1 1 0 0 0 0;
                                             1 1 0 0 1 1 0 0;
                                             1 0 1 0 1 0 1 0]
    #RM(2,3)
    @test parity_checks(ReedMuller(2,3)) == [1 1 1 1 1 1 1 1;
                                             1 1 1 1 0 0 0 0;
                                             1 1 0 0 1 1 0 0;
                                             1 0 1 0 1 0 1 0;
                                             1 1 0 0 0 0 0 0;
                                             1 0 1 0 0 0 0 0;
                                             1 0 0 0 1 0 0 0]
    #RM(3,3)
    @test parity_checks(ReedMuller(3,3)) == [1 1 1 1 1 1 1 1;
                                             1 1 1 1 0 0 0 0;
                                             1 1 0 0 1 1 0 0;
                                             1 0 1 0 1 0 1 0;
                                             1 1 0 0 0 0 0 0;
                                             1 0 1 0 0 0 0 0;
                                             1 0 0 0 1 0 0 0;
                                             1 0 0 0 0 0 0 0]
    #RM(2,4)
    @test parity_checks(ReedMuller(2,4)) == [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
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
