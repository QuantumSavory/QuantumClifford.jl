@testitem "Bravyi-Bacon-Shor and SHP Codes" tags=[:ecc] begin
    using QuantumClifford.ECC
    using QuantumClifford

    @testset "BBS Code from Hamming [7,4,3]" begin
        # As given in arXiv:2002.06257 Section 2.2
        A = [0 0 1 0 0 1 1;
             0 1 0 1 0 1 0;
             1 0 0 0 1 1 0;
             0 1 0 0 1 0 1;
             0 0 1 1 1 0 0;
             1 1 1 0 0 0 0;
             1 0 0 1 0 0 1]
        
        bbs = BravyiBaconShor(A)
        # N = 21, K = 4, D = 3
        @test code_n(bbs) == 21
        @test code_k(bbs) == 4
        
        # Table 1 implies 24 stabilizers (12 X-type, 12 Z-type)
        @test length(parity_checks(bbs)) == 24
        
        # Test CSS
        @test iscss(bbs)
        
        # Test exact number of stabilizers
        @test size(parity_matrix_x(bbs), 1) == 12
        @test size(parity_matrix_z(bbs), 1) == 12
    end
    
    @testset "SHP Code from Hamming [7,4,3]" begin
        # As given in arXiv:2002.06257 Section 3.3
        H = [1 1 0 1 1 0 0;
             1 0 1 1 0 1 0;
             0 1 1 1 0 0 1]
        
        shp = SubsystemHypergraphProduct(H, H)
        # N = 49, K = 16, D = 3
        @test code_n(shp) == 49
        @test code_k(shp) == 16
        
        # 24 stabilizer generators mentioned in the paper
        @test length(parity_checks(shp)) == 24
        
        # Test CSS
        @test iscss(shp)
        
        # Test sizes of X and Z parity check matrices
        @test size(parity_matrix_x(shp), 1) == 12
        @test size(parity_matrix_z(shp), 1) == 12
    end
end
