using WeakDepHelpers: @declare_method_is_in_extension, @declare_struct_is_in_extension

const heckeext_struct_docstring = "Implemented as a package extension with Hecke. Check the docs for the [Hecke extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Hecke.jl)"
const heckeext_function_docstring = heckeext_struct_docstring
const oscarext_struct_docstring = "Implemented as a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"
const oscarext_function_docstring = oscarext_struct_docstring

@declare_struct_is_in_extension QuantumClifford LPCode :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford LaCross :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford GeneralizedHyperGraphProduct :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford GeneralizedBicycle :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford ExtendedGeneralizedBicycle :QuantumCliffordHeckeExt (:Hecke,) heckeext_struct_docstring
@declare_struct_is_in_extension QuantumClifford DDimensionalSurface :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford DDimensionalToric :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford HomologicalProduct :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford DoubleHomologicalProduct :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford GeneralizedToric :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford TrivariateTricycle :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford BivariateBicycleViaPoly :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring
@declare_struct_is_in_extension QuantumClifford MultivariateMulticycle :QuantumCliffordOscarExt (:Oscar,) oscarext_struct_docstring

@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS two_block_group_algebra_code (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS generalized_bicycle_code_as_2bga (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS bicycle_code_as_2bga (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS Haah_cubic_code_as_2bga (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS honeycomb_color_code_as_2bga (:Hecke,) heckeext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS boundary_maps (:Oscar,) oscarext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS max_xy_exponents (:Oscar,) oscarext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS twobga_from_direct_product (:Oscar,) oscarext_function_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS twobga_from_fp_group (:Oscar,) oscarext_function_docstring
