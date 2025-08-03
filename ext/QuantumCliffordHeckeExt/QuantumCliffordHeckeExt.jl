module QuantumCliffordHeckeExt

using QECCore
import QECCore: code_n, code_s, code_k, rate, distance
using DocStringExtensions

import QuantumClifford
import QuantumClifford: gf2_row_echelon_with_pivots!
import LinearAlgebra
import Hecke: Group, GroupElem, AdditiveGroupElem,
    GroupAlgebra, GroupAlgebraElem, FqFieldElem, representation_matrix, dim, base_ring,
    multiplication_table, coefficients, abelian_group, group_algebra, rand, gens, order,
    is_commutative, FqPolyRingElem, residue_ring, coeff, zero_matrix, mod1, lift, ZZ, gen,
    matrix, ncols, nrows, degree, gcd, polynomial_ring
import Nemo
import Nemo: characteristic, matrix_repr, GF, ZZ, lift

import QuantumClifford.ECC: iscss, parity_checks,
    two_block_group_algebra_codes, generalized_bicycle_codes, bicycle_codes, check_repr_commutation_relation,
    haah_cubic_codes, honeycomb_color_codes, check_repr_regular_linear

import QECCore: AbstractECC, CSS,
    hgp, code_k, code_n, code_s, parity_matrix_x, parity_matrix_z, parity_matrix_xz

# exported from extension so that Documenter.jl sees them when autogenerating API lists
export hgp, two_block_group_algebra_codes, generalized_bicycle_codes, bicycle_codes, haah_cubic_codes,
    LPCode, LiftedCode, honeycomb_color_codes, LaCross, GeneralizedBicycleCode

include("util.jl")
include("types.jl")
include("lifted.jl")
include("lacross.jl")
include("lifted_product.jl")
include("generalized_bicycle.jl")

end # module
