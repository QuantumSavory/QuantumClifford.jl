module QuantumCliffordHeckeExt

using QECCore
import QECCore: code_n, code_s, code_k, rate, distance
using DocStringExtensions

import QuantumClifford, LinearAlgebra
import Hecke: Group, GroupElem, AdditiveGroupElem,
    GroupAlgebra, GroupAlgebraElem, FqFieldElem, representation_matrix, dim, base_ring,
    multiplication_table, coefficients, abelian_group, group_algebra, rand, is_commutative
import Nemo
import Nemo: characteristic, matrix_repr, GF, ZZ, lift

import QuantumClifford.ECC: AbstractECC, CSS,
    hgp, code_k, code_n, code_s, iscss, parity_checks, parity_checks_x, parity_checks_z, parity_checks_xz,
    two_block_group_algebra_codes, generalized_bicycle_codes, bicycle_codes, check_repr_commutation_relation,
    haah_cubic_codes

include("util.jl")
include("types.jl")
include("lifted.jl")
include("lifted_product.jl")

end # module
