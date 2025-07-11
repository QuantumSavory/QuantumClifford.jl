module QECCore

using SparseArrays
using LinearAlgebra
using Combinatorics
using Graphs

using DocStringExtensions

using DocStringExtensions

# interfaces
export distance, parity_matrix, code_n, code_s, code_k, parity_matrix_x, parity_matrix_z, parity_matrix_xz, rate, metacheck_matrix_x, metacheck_matrix_z, metacheck_matrix
export AbstractECC, AbstractQECC, AbstractCECC, AbstractCSSCode, AbstractDistanceAlg

# QEC Codes
export Perfect5, Cleve8, Gottesman

# CSS Codes

export Toric, Bitflip3, Phaseflip3, Shor9, Steane7, Surface, CSS, QuantumReedMuller, Triangular488, Triangular666, QuantumTannerGraphProduct, CyclicQuantumTannerGraphProduct

# Classical Codes
export RepCode, ReedMuller, RecursiveReedMuller, Golay, Hamming

include("interface.jl")
include("codes/util.jl")

# Classical Codes
include("codes/classical/hamming.jl")
include("codes/classical/reptetion.jl")
include("codes/classical/golay.jl")

# Quantum Codes
include("codes/quantum/css.jl")
include("codes/quantum/fivequbit.jl")
include("codes/quantum/toric.jl")
include("codes/quantum/clevecode.jl")
include("codes/quantum/shorcode.jl")
include("codes/quantum/steanecode.jl")
include("codes/quantum/surface.jl")
include("codes/quantum/bitflipcode.jl")
include("codes/quantum/gottesman.jl")
include("codes/quantum/color_codes.jl")
include("codes/quantum/quantumtannergraphproduct.jl")

# Reed-Muller Codes
include("codes/classical/reedmuller.jl")
include("codes/classical/recursivereedmuller.jl")
include("codes/quantum/quantumreedmuller.jl")

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
