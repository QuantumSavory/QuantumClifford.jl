module QECCoreNemoExt

using QECCore
using DocStringExtensions

import Nemo
import Nemo: GF, matrix, rank, transpose, finite_field, GF, polynomial_ring, evaluate,
    FqFieldElem, FqPolyRingElem, degree, is_irreducible, gcd, derivative, matrix, inv,
    is_zero, coeff, is_monic, one

import QECCore: code_k, parity_matrix_x, parity_matrix_z, parity_matrix, generator_polynomial

import Random
import Random: MersenneTwister, GLOBAL_RNG, AbstractRNG, rand

import QECCore: random_Goppa_code, code_k, code_n

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

include("goppa.jl")

end # module
