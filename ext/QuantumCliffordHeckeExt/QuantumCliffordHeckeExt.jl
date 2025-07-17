module QuantumCliffordHeckeExt

using QECCore
import QECCore: code_n, code_s, code_k, rate, distance
using DocStringExtensions

import QuantumClifford
import LinearAlgebra
import LinearAlgebra: rank
import Hecke: Group, GroupElem, AdditiveGroupElem,
    GroupAlgebra, GroupAlgebraElem, FqFieldElem, representation_matrix, dim, base_ring,
    multiplication_table, coefficients, abelian_group, group_algebra, rand, gens, order,
    is_commutative, group
import Nemo
import Nemo: characteristic, matrix_repr, GF, ZZ, lift, matrix

import QuantumClifford.ECC: iscss, parity_checks, parity_matrix,
    two_block_group_algebra_codes, generalized_bicycle_codes, bicycle_codes, check_repr_commutation_relation,
    haah_cubic_codes, honeycomb_color_codes, check_repr_regular_linear

import QECCore: AbstractECC, CSS,
    hgp, code_k, code_n, code_s, parity_matrix_x, parity_matrix_z, parity_matrix_xz, parity_matrix

# exported from extension so that Documenter.jl sees them when autogenerating API lists
export hgp, two_block_group_algebra_codes, generalized_bicycle_codes, bicycle_codes, haah_cubic_codes,
    LPCode, LiftedCode, honeycomb_color_codes

include("util.jl")
include("types.jl")
include("lifted.jl")
include("lifted_product.jl")

end # module
