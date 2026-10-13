module QECCoreNemoExt

using QECCore
using DocStringExtensions

import Nemo
import Nemo: GF, gen, matrix, rank, transpose, polynomial_ring, evaluate, FqFieldElem,
    FqPolyRingElem, degree, is_irreducible, gcd, derivative, inv, coeff, is_monic, one,
    minpoly, lcm, is_zero, finite_field
import QECCore: AbstractPolynomialCode, BCH, code_k, code_n, generator_polynomial,
    parity_matrix, parity_matrix_x, parity_matrix_z, random_Goppa_code

import Random
import Random: MersenneTwister, GLOBAL_RNG, AbstractRNG, rand

function QECCore.code_k(c::AbstractCSSCode)
    n = code_n(c)
    F₂ = GF(2)
    rank_Hx = rank(matrix(F₂, parity_matrix_x(c)))
    rank_Hz = rank(matrix(F₂, parity_matrix_z(c)))
    return n - rank_Hx - rank_Hz
end

function QECCore.code_k(c::AbstractCECC)
    H = Matrix(parity_matrix(c))
    n = code_n(c)
    F₂ = GF(2)
    rank_H = rank(matrix(F₂, H))
    return n - rank_H
end

include("bch.jl")
include("goppa.jl")

end # module
