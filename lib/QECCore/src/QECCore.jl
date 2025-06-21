module QECCore

using SparseArrays
using LinearAlgebra
using Combinatorics

# interfaces
export distance, parity_matrix, code_n, code_s, code_k, parity_matrix_x, parity_matrix_z, rate
export AbstractECC, AbstractQECC, AbstractCECC, AbstractCSSCode, AbstractDistanceAlg

# QEC Codes
export Perfect5, Cleve8

# CSS Codes
export Toric, Bitflip3, Phaseflip3, Shor9, Steane7, Surface, CSS

# Classical Codes
export RepCode

include("interface.jl")
include("codes/util.jl")
include("codes/css.jl")
include("codes/fivequbit.jl")
include("codes/reptetion.jl")
include("codes/toric.jl")
include("codes/clevecode.jl")
include("codes/shorcode.jl")
include("codes/steanecode.jl")
include("codes/surface.jl")
include("codes/bitflipcode.jl")

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f == parity_matrix_x || exc.f == parity_matrix_z
                print(io, "\n Hint: This error usually occurs when trying to call `parity_matrix_x` or `parity_matrix_z` on a non-CSS quantum code.")
            end
        end
    end
end

end
