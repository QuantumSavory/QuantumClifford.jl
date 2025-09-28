module QuantumCliffordOscarExt

using DocStringExtensions

import LinearAlgebra
import LinearAlgebra: rank
import QuantumClifford
import QuantumClifford: Stabilizer
import Nemo
import Nemo: FqFieldElem
import Hecke: group_algebra, GF, abelian_group, gens, quo, one, GroupAlgebra,
   GroupAlgebraElem, direct_product, sub, ZZ, lift, polynomial_ring
import Oscar
import Oscar: free_group, small_group_identification, describe, order, FPGroupElem, FPGroup,
    BasicGAPGroupElem, DirectProductGroup, cyclic_group, free_module, hom, transpose, tensor_product,
    chain_complex, total_complex, map, summands, MatElem, matrix, nrows, ncols, kernel, dim, range, image,
    base_ring, ComplexOfMorphisms, coefficients, zero_matrix, hcat, circshift, size, zeros, enumerate,
    kronecker_product, FqMatrix, identity_matrix, iszero, FqPolyRingElem, laurent_polynomial_ring,
    hnf_with_transform, ideal, intersect, ==, is_coprime, quo, groebner_basis, length, FqMPolyRingElem,
    first, length, MPolyQuoRingElem, FqMPolyRingElem, modulus, ideal, monomials, terms, coeff, degree, mod,
    koszul_matrix
import Oscar.Generic.MatSpaceElem
import Oscar.Generic.DirectSumModule
import Oscar.Generic.LaurentMPolyWrap
import Oscar.Generic.exponent_vectors
import Oscar.IdealGens
import Combinatorics: combinations

import QuantumClifford.ECC: two_block_group_algebra_codes, twobga_from_direct_product, twobga_from_fp_group,
    boundary_maps, max_xy_exponents

import QECCore: AbstractECC, CSS, RepCode, AbstractCSSCode,
    hgp, code_k, code_n, code_s, distance, parity_matrix_x, parity_matrix_z, parity_matrix_xz, parity_matrix,
    metacheck_matrix_x, metacheck_matrix_z, metacheck_matrix

# exported from extension so that Documenter.jl sees them when autogenerating API lists
export twobga_from_direct_product, twobga_from_fp_group, DDimensionalSurfaceCode, DDimensionalToricCode, boundary_maps,
    HomologicalProductCode, DoubleHomologicalProductCode, GeneralizedToricCode, TrivariateTricycleCode

include("types.jl")
include("direct_product.jl")
include("generalized_toric.jl")
include("group_presentation.jl")
include("d_dimensional_codes.jl")
include("trivariate_tricycle.jl")
include("multivariate_multicycle.jl")
include("homological_product_codes.jl")

end # module
