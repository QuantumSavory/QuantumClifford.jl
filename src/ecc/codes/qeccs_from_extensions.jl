using WeakDepHelpers: @declare_method_is_in_extension, @declare_struct_is_in_extension

const heckeext_struct_docstring = "Implemented as a package extension with Hecke. Check the docs for the [Hecke extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Hecke.jl)"
const heckeext_function_docstring = heckeext_struct_docstring
const oscarext_struct_docstring = "Implemented as a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"
const oscarext_function_docstring = oscarext_struct_docstring

@declare_struct_is_in_extension QuantumClifford LPCode :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford LaCross :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford GeneralizedHyperGraphProductCode :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford GeneralizedBicycleCode :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford ExtendedGeneralizedBicycleCode :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford DDimensionalSurfaceCode :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford DDimensionalToricCode :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford HomologicalProductCode :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford DoubleHomologicalProductCode :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford GeneralizedToricCode :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford TrivariateTricycleCode :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford BivariateBicycleCode :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring

@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS two_block_group_algebra_codes (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS generalized_bicycle_codes (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS bicycle_codes (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS haah_cubic_codes (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS honeycomb_color_codes (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS boundary_maps (:Oscar,) oscarext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS max_xy_exponents (:Oscar,) oscarext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS twobga_from_direct_product (:Oscar,) oscarext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS twobga_from_fp_group (:Oscar,) oscarext_function_docstring
