module QuantumCliffordHeckeExt

using DocStringExtensions

import QuantumClifford, LinearAlgebra
import Hecke: Group, GroupElem, AdditiveGroup, AdditiveGroupElem,
    GroupAlgebra, GroupAlgebraElem, FqFieldElem, representation_matrix, dim, base_ring,
    multiplication_table, coefficients, abelian_group, group_algebra
import Nemo: characteristic, matrix_repr, GF, ZZ

import QuantumClifford.ECC: AbstractECC, CSS, ClassicalCode,
    hgp, code_k, code_n, code_s, iscss, parity_checks, parity_checks_x, parity_checks_z, parity_checks_xz,
    two_block_group_algebra_codes, generalized_bicycle_codes, bicycle_codes

include("types.jl")
include("lifted.jl")
include("lifted_product.jl")

end # module
