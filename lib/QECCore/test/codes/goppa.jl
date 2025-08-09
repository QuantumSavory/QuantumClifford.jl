@testitem "ECC Golay" begin
    using QECCore
    using LinearAlgebra
    using Nemo: finite_field, GF, polynomial_ring, evaluate, FqFieldElem, degree, is_irreducible, gcd, derivative, matrix, inv
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC

    @testset "Testing Goppa codes properties" begin
        for m in 3:12
            for t in 2:m
                gx = generator_polynomial(GoppaCode(m, t))
                # Goppa generator polynomial, `g(x)` is an irreducible polynomial.
                @test is_irreducible(gx) == true
                # The Goppa polynomial, `g(x)`, is square free, thus no repeated roots, if the gcd with the derivative is 1.
                @test gcd(derivative(gx), gx) == 1
                mat = matrix(GF(2), parity_matrix(GoppaCode(m, t)))
                computed_rank = rank(mat)
                # For parity check matrix matrix `H = XY` over `GF(2ᵐ)`, maximum of `m * t` rows are linearly independent, hence `Rank(XY) ≤ m * t`[singh2019code](@cite).
                @test computed_rank <= m * t
            end
        end
    end
end
