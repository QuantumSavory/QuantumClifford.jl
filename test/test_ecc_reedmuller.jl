using Test
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC

@testset "Testing common examples RM(r,m) codes [raaphorst2003reed](@cite) and [abbe2020reed](cite)" begin
    
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