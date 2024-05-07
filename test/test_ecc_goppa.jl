using Test
using LinearAlgebra
using Nemo: finite_field, GF, polynomial_ring, evaluate, FqFieldElem, degree, is_irreducible, gcd, derivative, matrix, inv
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, generator_polynomial, Goppa

function designed_distance(matrix, t)
    for row in eachrow(matrix)
        count = sum(row)
        # The minimum Hamming distance of Goppa code, `d(Γ(L, g)) ≥ t + 1` [singh2019code](@cite).
        if count >= t + 1
            return true
        end
    end
    return false
end

@testset "Testing Goppa codes properties" begin
    m_cases = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    for m in m_cases
        for t in 2:m
            gx = generator_polynomial(Goppa(m, t))
            # Goppa generator polynomial, `g(x)` is an irreducible polynomial.
            @test is_irreducible(gx) == true
            # The Goppa polynomial, `g(x)`, is square free, thus no repeated roots, if the gcd with the derivative is 1.
            @test gcd(derivative(gx), gx) == 1
            @test designed_distance(parity_checks(Goppa(m, t)), t) == true
            mat = matrix(GF(2), parity_checks(Goppa(m, t)))
            computed_rank = rank(mat)
            # For parity check matrix matrix `H = XY` over `GF(2ᵐ)`, maximum of `m * t` rows are linearly independent, hence `Rank(XY) ≤ m * t`[singh2019code](@cite).
            @test computed_rank <= m * t
        end
    end
end
