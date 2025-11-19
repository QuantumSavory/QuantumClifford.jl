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
    matrix, ncols, nrows, degree, EuclideanRingResidueRingElem, quo, parent, zero, gcd,
    polynomial_ring, characteristic, isone, mod, factor, zeros
import Hecke.Generic.MatSpaceElem
import Nemo
import Nemo: characteristic, matrix_repr, GF, ZZ, lift, rank

import QuantumClifford.ECC: iscss, parity_checks,
    two_block_group_algebra_codes, generalized_bicycle_codes_as_2bga, bicycle_codes_as_2bga, check_repr_commutation_relation,
    Haah_cubic_codes_as_2bga, honeycomb_color_codes_as_2bga, check_repr_regular_linear, random_qc_ghp_code_matrix_A

import QECCore: AbstractQECC, CSS, AbstractCSSCode,
    hgp, code_k, code_n, code_s, parity_matrix_x, parity_matrix_z, parity_matrix_xz

import Random
import Random: AbstractRNG, default_rng, randperm

# exported from extension so that Documenter.jl sees them when autogenerating API lists
export hgp, two_block_group_algebra_codes, generalized_bicycle_codes_as_2bga, bicycle_codes_as_2bga, Haah_cubic_codes_as_2bga,
    LPCode, LiftedCode, honeycomb_color_codes_as_2bga, LaCross, GeneralizedBicycleCode, ExtendedGeneralizedBicycleCode,
    GeneralizedHyperGraphProductCode

include("util.jl")
include("types.jl")
include("lifted.jl")
include("lacross.jl")
include("lifted_product.jl")
include("generalized_bicycle.jl")
include("extended_generalized_bicycle.jl")
include("generalized_hypergraph_product.jl")

end # module
