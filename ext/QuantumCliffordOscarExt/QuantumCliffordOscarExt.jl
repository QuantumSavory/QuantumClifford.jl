module QuantumCliffordOscarExt

using DocStringExtensions

import LinearAlgebra
import LinearAlgebra: rank
import Nemo
import Nemo: FqFieldElem
import Hecke: group_algebra, GF, abelian_group, gens, quo, one, GroupAlgebra,
   GroupAlgebraElem, direct_product, sub, ZZ, lift
import Oscar
import Oscar: free_group, small_group_identification, describe, order, FPGroupElem, FPGroup,
    BasicGAPGroupElem, DirectProductGroup, cyclic_group, free_module, hom, transpose, tensor_product,
    chain_complex, total_complex, map, summands, MatElem, matrix, nrows, ncols
import Oscar.Generic.DirectSumModule

import QuantumClifford.ECC: two_block_group_algebra_codes, twobga_from_direct_product, twobga_from_fp_group, d_dimensional_surface_codes

import QECCore: AbstractECC, CSS, RepCode,
    hgp, code_k, code_n, code_s, parity_matrix_x, parity_matrix_z, parity_matrix_xz, parity_matrix

# exported from extension so that Documenter.jl sees them when autogenerating API lists
export twobga_from_direct_product, twobga_from_fp_group, d_dimensional_surface_codes

include("types.jl")
include("direct_product.jl")
include("group_presentation.jl")
include("d_dimensional_surface_codes.jl")

end # module
