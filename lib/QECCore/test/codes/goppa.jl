@testitem "ECC Goppa" begin
    using Random
    using QECCore
    using LinearAlgebra
    import Nemo
    using Nemo: finite_field, GF, polynomial_ring, evaluate, FqFieldElem, degree, is_irreducible, gcd, derivative, matrix, inv, is_monic, echelon_form, transpose, nullspace, iszero
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC
    using Random: MersenneTwister, GLOBAL_RNG, AbstractRNG, rand
    using QECCore: random_Goppa_code

    @testset "Testing Goppa codes properties" begin
        rng = MersenneTwister(123)
        for m in 3:12
            for t in 2:m-1
                ga = random_Goppa_code(rng, m, t)
                gx = ga.g
                # Goppa generator polynomial, `g(x)` is an irreducible polynomial.
                @test is_irreducible(gx) && is_monic(gx)
                # The Goppa polynomial, `g(x)`, is square free, thus no repeated roots, if the gcd with the derivative is 1.
                @test gcd(derivative(gx), gx) == 1
                mat = matrix(GF(2), parity_matrix(ga))
                @test size(mat, 1) <= size(mat, 2)
                computed_rank = rank(mat)
                # For parity check matrix matrix `H = XY` over `GF(2ᵐ)`, maximum of `m * t` rows are linearly independent, hence `Rank(XY) ≤ m * t`[singh2019code](@cite).
                @test computed_rank <= m * t
                n = length(ga.L)
                @test  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) >= n - m*t
                @test degree(ga.g) ≤ length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) # degenerate cases included.
            end
        end
    end

    # [12, 4, ≥ 5]
    m = 4
    t = 2
    F, α = finite_field(2, m, :α)
    R, x = polynomial_ring(F, :x)
    g = x^2 + x + α^3
    L = [α^i for i in 2:13]
    ga = GoppaCode(m, t, g, L)
    n = length(L)
    k = n - m*t 
    d = 2t + 1 
    H = parity_matrix(ga)
    @test  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) >= n - m*t
    # From https://surface.syr.edu/cgi/viewcontent.cgi?article=1846&context=honors_capstone
    ref_mat = [0 0 1 0 1 0 0 0 0 0 0 1;
               0 1 1 1 0 0 0 0 1 0 0 0;
               0 0 0 0 1 1 0 1 0 1 1 1;
               1 1 0 0 0 1 1 1 0 0 0 0;
               1 1 1 1 0 0 1 0 1 1 1 1;
               0 0 1 0 0 0 0 0 1 1 1 0;
               0 1 0 0 0 0 0 1 0 1 0 1;
               0 1 1 1 1 0 1 1 1 0 0 1];
    @test echelon_form(matrix(GF(2), H)) == echelon_form(matrix(GF(2), ref_mat))
    @test gcd(derivative(ga.g), ga.g) == 1
    mat = matrix(GF(2), parity_matrix(ga))
    computed_rank = rank(mat)
    @test computed_rank <= m * t
    n = length(ga.L)
    @test  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) >= n - m*t

    # From https://doc.sagemath.org/html/en/reference/coding/sage/coding/goppa_code.html
    # [55, 16, d]
    m = 6
    F, α = finite_field(2, m, :α)
    R, x = polynomial_ring(F, :x)
    t = 9
    g = x^t + 1
    ga = GoppaCode(m, t, g)
    @test length(ga.L) == 55 &&  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) == 16
    @test gcd(derivative(ga.g), ga.g) == 1
    mat = matrix(GF(2), parity_matrix(ga))
    computed_rank = rank(mat)
    @test computed_rank <= m * t
    n = length(ga.L)
    @test  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) >= n - m*t
    # [8, 2, d]
    m = 3
    F, α = finite_field(2, m, :α)
    R, x = polynomial_ring(F, :x)
    t = 2
    g = x^t + x + 1
    ga = GoppaCode(m, t, g)
    @test length(ga.L) == 8 &&  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) == 2
    @test gcd(derivative(ga.g), ga.g) == 1
    mat = matrix(GF(2), parity_matrix(ga))
    computed_rank = rank(mat)
    @test computed_rank <= m * t
    n = length(ga.L)
    @test  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) >= n - m*t

    # From example 1.9 of https://arxiv.org/pdf/1907.12754
    m = 4
    t = 2
    F, α = finite_field(2, m, :α)
    R, z = polynomial_ring(F, :z)
    g = z^2 + α^7*z + 1
    L = [α^i for i in 2:13]
    ga = GoppaCode(m, t, g, L)
    n = length(L)
    k = n - m*t 
    d = 2t + 1 
    H = parity_matrix(ga)
    @test  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) >= n - m*t
    mat = matrix(GF(2), parity_matrix(ga))
    computed_rank = rank(mat)
    @test computed_rank <= m * t
    ref_mat = [0 1 0 1 0 0 1 1 0 1 1 0;
               1 1 1 0 0 0 1 0 0 1 0 0;
               0 1 0 0 1 0 1 1 1 0 0 1;
               1 0 1 1 1 0 0 0 0 1 1 1;
               1 0 0 0 0 1 1 1 0 1 1 0;
               1 0 0 1 1 0 0 0 1 0 1 1;
               1 1 1 0 1 0 0 1 1 0 1 0;
               1 1 1 0 1 0 1 0 1 1 1 0]
    @test echelon_form(matrix(GF(2), H)) == echelon_form(matrix(GF(2), ref_mat))
    ref_G = [1 1 1 1 0 1 0 1 0 1 0 0;
             0 1 0 0 1 1 1 1 0 0 1 0;
             0 0 1 0 1 0 1 1 1 0 0 0;
             0 1 0 1 0 0 1 1 0 0 0 1];
    G = transpose(nullspace(mat)[2])
    @test iszero(G*transpose(mat))
    @test echelon_form(G) == echelon_form(matrix(GF(2), ref_G))

    # From https://crypto-kantiana.com/elena.kirshanova/talks/Talk_McEliece.pdf
    t = 2
    m = 3
    F, α = finite_field(2, m, :α)
    R, x = polynomial_ring(F, :x)
    g = x^t + x + 1
    L = [F(0), F(1), α, α^2, α + 1, α^2 + α, α^2 + α + 1, α^2 + 1]
    ga = GoppaCode(m, t, g, L)
    H = parity_matrix(ga)
    ref_H = [1 1 0 0 0 0 0 0;
             0 0 0 1 0 1 1 1;
             0 0 1 1 1 0 0 1;
             0 1 1 1 1 1 1 1;
             0 0 1 0 1 1 0 1;
             0 0 0 1 1 1 1 0];
    mat = matrix(GF(2), parity_matrix(ga))
    computed_rank = rank(mat)
    @test computed_rank <= m * t
    n = length(ga.L)
    @test  length(ga.L) - rank(matrix(GF(2), parity_matrix(ga))) >= n - m*t
    @test echelon_form(matrix(GF(2), H)) == echelon_form(matrix(GF(2), ref_H))
end
