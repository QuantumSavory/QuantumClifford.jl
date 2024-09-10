module QuantumCliffordHeckeExt

using DocStringExtensions

import QuantumClifford, LinearAlgebra
import AbstractAlgebra: Group, GroupElem, AdditiveGroup, AdditiveGroupElem
import Hecke: GroupAlgebra, GroupAlgebraElem, FqFieldElem, representation_matrix, dim, base_ring,
    multiplication_table, coefficients, abelian_group, group_algebra
import Base: adjoint
import Nemo: characteristic, lift, matrix_repr, GF, ZZ

import QuantumClifford.ECC: AbstractECC, CSS, ClassicalCode,
    hgp, code_k, code_n, code_s, iscss, parity_checks, parity_checks_x, parity_checks_z, parity_checks_xz

include("types.jl")
include("lifted.jl")
include("lifted_product.jl")

end # module