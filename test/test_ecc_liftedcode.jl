@testitem "ECC LiftedCode" tags=[:ecc, :ecc_bespoke] begin
    using Hecke
    using Hecke: group_algebra, GF, abelian_group, gens, one, representation_matrix
    using QuantumClifford.ECC: LiftedCode, code_k, code_n, code_s, parity_matrix

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
