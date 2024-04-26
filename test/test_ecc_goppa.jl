using Test
using Nemo: finite_field, GF, polynomial_ring, evaluate, FqFieldElem, degree, is_irreducible
using AbstractAlgebra
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, generator_polynomial

function designed_distance(matrix, t)
    for row in eachrow(matrix)
        count = sum(row)
        if count >= t || count >= t - 1 || count >= t - 2 || count >= t - 3
            return true
        end
    end
    return false
end

@testset "Testing Goppa codes properties" begin
    n_cases = [8, 16, 32, 64, 128, 256]
    for n in n_cases
        for t in 3:7
            r = ceil(Int, log2(n))
            GF2ͬ, o = finite_field(2, r, "o")
            k = GF(2, r)
            po, b = polynomial_ring(k)
            gx = generator_polynomial(Goppa(n, t))
            L = FqFieldElem[]
            i = 0 
            while length(L) != n
                if evaluate(gx, o ^ i) != 0
                    L = [L; evaluate(gx, o^i)]
                end
                i += 1
            end
            @test is_irreducible(gx) == true
            @test degree(gx) == t
            @test gcd(b - L[t], evaluate(gx, b)) == 1
            @test designed_distance(parity_checks(Goppa(n, t)), t) == true
        end
    end
end
