module QuantumCliffordOscarExt

using DocStringExtensions

import Nemo
import Nemo: FqFieldElem
import Hecke: group_algebra, GF, abelian_group, gens, quo, one, GroupAlgebra,
   GroupAlgebraElem, direct_product, sub
import Oscar
import Oscar: free_group, small_group_identification, describe, order, FPGroupElem, FPGroup,
    BasicGAPGroupElem, DirectProductGroup, cyclic_group

import QuantumClifford.ECC: two_block_group_algebra_codes

include("types.jl")
include("direct_product.jl")
include("group_presentation.jl")

end # module
