using QuantumClifford: register_method_error_hint, WeakDepMissingError

const heckeext_struct_docstring = "Implemented as a package extension with Hecke. Check the docs for the [Hecke extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Hecke.jl)"
const heckeext_function_docstring = heckeext_struct_docstring
const oscarext_struct_docstring = "Implemented as a package extension with Oscar. Check the docs for the [Oscar extension](http://qc.quantumsavory.org/stable/ECC_API/#Implemented-in-an-extension-requiring-Oscar.jl)"
const oscarext_function_docstring = oscarext_struct_docstring

const structs_implemented_in_extensions = [
    (:LPCode, :QuantumCliffordHeckeExt, heckeext_struct_docstring, (:Hecke,)),
    (:LaCross, :QuantumCliffordHeckeExt, heckeext_struct_docstring, (:Hecke,)),
    (:GeneralizedHyperGraphProductCode, :QuantumCliffordHeckeExt, heckeext_struct_docstring, (:Hecke,)),
    (:GeneralizedBicycle, :QuantumCliffordHeckeExt, heckeext_struct_docstring, (:Hecke,)),
    (:ExtendedGeneralizedBicycle, :QuantumCliffordHeckeExt, heckeext_struct_docstring, (:Hecke,)),
    (:DDimensionalSurface, :QuantumCliffordOscarExt, oscarext_struct_docstring, (:Oscar,)),
    (:DDimensionalToric, :QuantumCliffordOscarExt, oscarext_struct_docstring, (:Oscar,)),
    (:HomologicalProduct, :QuantumCliffordOscarExt, oscarext_struct_docstring, (:Oscar,)),
    (:DoubleHomologicalProduct, :QuantumCliffordOscarExt, oscarext_struct_docstring, (:Oscar,)),
    (:GeneralizedToric, :QuantumCliffordOscarExt, oscarext_struct_docstring, (:Oscar,)),
    (:TrivariateTricycle, :QuantumCliffordOscarExt, oscarext_struct_docstring, (:Oscar,)),
    (:BivariateBicycleCodeViaPoly, :QuantumCliffordOscarExt, oscarext_struct_docstring, (:Oscar,)),
    (:MultivariateMulticycle, :QuantumCliffordOscarExt, oscarext_struct_docstring, (:Oscar,)),
]

for (struct_name, extension_name, struct_docstring, deps) in structs_implemented_in_extensions
    expr = quote
        function $struct_name(args...; kwargs...)
            ext = Base.get_extension(QuantumClifford, $(Base.Meta.quot(extension_name)))
            if isnothing(ext)
                throw(WeakDepMissingError($(Base.Meta.quot(struct_name)), $(deps)))
            end
            return ext.$struct_name(args...; kwargs...)
        end
        @doc $(struct_docstring) $(struct_name)
    end
    eval(expr)
end

const functions_implemented_in_extensions = [
    (:two_block_group_algebra_code, :QuantumCliffordHeckeExt, heckeext_function_docstring, (:Hecke,)),
    (:generalized_bicycle_code_as_2bga, :QuantumCliffordHeckeExt, heckeext_function_docstring, (:Hecke,)),
    (:bicycle_code_as_2bga, :QuantumCliffordHeckeExt, heckeext_function_docstring, (:Hecke,)),
    (:Haah_cubic_code_as_2bga, :QuantumCliffordHeckeExt, heckeext_function_docstring, (:Hecke,)),
    (:honeycomb_color_code_as_2bga, :QuantumCliffordHeckeExt, heckeext_function_docstring, (:Hecke,)),
    (:boundary_maps, :QuantumCliffordOscarExt, oscarext_function_docstring, (:Oscar,)),
    (:max_xy_exponents, :QuantumCliffordOscarExt, oscarext_function_docstring, (:Oscar,)),
    (:twobga_from_direct_product, :QuantumCliffordOscarExt, oscarext_function_docstring, (:Oscar,)),
    (:twobga_from_fp_group, :QuantumCliffordOscarExt, oscarext_function_docstring, (:Oscar,)),
]

for (function_name, extension_name, function_docstring, deps) in functions_implemented_in_extensions
    expr = quote
        function $function_name end
    end
    fe = eval(expr)
    eval(:(@doc $function_docstring $function_name))
    register_method_error_hint(fe, deps)
end
