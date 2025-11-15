const BIG_INT_MINUS_ONE = Ref{BigInt}()
const BIG_INT_TWO = Ref{BigInt}()
const BIG_INT_FOUR = Ref{BigInt}()

struct WeakDepMissingError <: Exception
    name::Symbol
    deps::Tuple{Vararg{Symbol}}
end

function Base.showerror(io::IO, e::WeakDepMissingError)
    hl_name = highlight(string(e.name))
    hl_deps = highlight(string(e.deps...))
    hl_import_deps = highlight(string("import ", join(e.deps, ", ")))
    print(io, styled"{info:`$(hl_name)` depends on the package(s) `$(hl_deps)` but you have not installed or imported them yet. Immediately after an `$(hl_import_deps)`, `$(hl_name)` will be available.}")
end

const METHOD_ERROR_HINTS_FOR_EXTENSIONS = Dict{Function, Tuple{Vararg{Symbol}}}()

function register_method_error_hint(f::Function, deps::Tuple{Vararg{Symbol}})
    METHOD_ERROR_HINTS_FOR_EXTENSIONS[f] = deps
end

function method_error_hint_callback(io, exc, argtypes, kwargs)
    if haskey(METHOD_ERROR_HINTS_FOR_EXTENSIONS, exc.f)
        deps = METHOD_ERROR_HINTS_FOR_EXTENSIONS[exc.f]
        print(io, styled"\n{bold:{info:HINT: }}")
        showerror(io, WeakDepMissingError(nameof(exc.f), deps))
    end
end

function __init__()
    BIG_INT_MINUS_ONE[] = BigInt(-1)
    BIG_INT_TWO[] = BigInt(2)
    BIG_INT_FOUR[] = BigInt(4)

    # Register error hint for the `project!` method for GeneralizedStabilizer
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            method_error_hint_callback(io, exc, argtypes, kwargs)
            if exc.f === project! && argtypes[1] <: GeneralizedStabilizer
                print(io, """
                \nThe method `project!` is not appropriate for use with`GeneralizedStabilizer`.
                You probably are looking for `projectrand!`.
                `project!` in this library is a low-level "linear algebra" method to verify
                whether a measurement operator commutes with a set of stabilizers, and to
                potentially simplify the tableau and provide the index of the anticommuting
                term in that tableau. This linear algebra operation is not defined for
                `GeneralStabilizer` as there is no single tableau to provide an index into.""")
            elseif exc.f === ECC.distance && length(argtypes)==1
                print(io,"""
                \nThe distance for this code is not in our database. Consider using the MIP-based method:
                `import JuMP, HiGHS; distance(code, DistanceMIPAlgorithm(solver=HiGHS))` or another MIP solver""")
            elseif exc.f === ECC.distance && length(argtypes)==2 && argtypes[2]===ECC.DistanceMIPAlgorithm
                print(io,"""\nPlease first import `JuMP` to make MIP-based distance calculation available.""")
            end
        end
    end
end
