@testitem "ECC throws" tags=[:ecc, :ecc_base] begin

    using Hecke
    using Hecke: group_algebra, GF, abelian_group, gens, one, representation_matrix
    using QuantumClifford.ECC: ReedMuller, BCH, RecursiveReedMuller, Golay, Triangular488, Triangular666, Hamming, LiftedCode, code_k, code_n, code_s, parity_matrix

    @test_throws ArgumentError ReedMuller(-1, 3)
    @test_throws ArgumentError ReedMuller(1, 0)
    @test_throws ArgumentError ReedMuller(4, 2)

    @test_throws ArgumentError BCH(2, 2)
    @test_throws ArgumentError BCH(-2, 2)
    @test_throws ArgumentError BCH(2, -3)
    @test_throws ArgumentError BCH(3, 3)
    @test_throws ArgumentError BCH(3, 4)
    @test_throws ArgumentError BCH(4, 4)
    @test_throws ArgumentError BCH(4, 100)

    @test_throws ArgumentError RecursiveReedMuller(-1, 3)
    @test_throws ArgumentError RecursiveReedMuller(1, 0)
    @test_throws ArgumentError RecursiveReedMuller(4, 2)

    @test_throws ArgumentError Golay(21)

    @test_throws ArgumentError Hamming(1)

    @test_throws ArgumentError Triangular488(8)
    @test_throws ArgumentError Triangular488(1)
    @test_throws ArgumentError Triangular666(8)
    @test_throws ArgumentError Triangular666(1)

    l = 12; GA = group_algebra(GF(2), abelian_group(l)); x = gens(GA)[]
    B = reshape([1 + x + x^3 + x^6], (1, 1))
    c = LiftedCode(B, repr = representation_matrix)
    @test_nowarn parity_matrix(c)
    @test_nowarn code_n(c)
    @test_nowarn code_k(c)
    @test_nowarn code_s(c)
    base_matrix = [0 0 0 0; 0 1 2 5; 0 6 3 1]; l = 3;
    c = LiftedCode(base_matrix, l);
    @test_nowarn parity_matrix(c)
    @test_nowarn code_n(c)
    @test_nowarn code_k(c)
    @test_nowarn code_s(c)
end
