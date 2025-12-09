module QECCore

using SparseArrays
using LinearAlgebra
using Combinatorics
using Graphs
using Random: GLOBAL_RNG, AbstractRNG, randperm, rand, MersenneTwister, randperm

using DocStringExtensions

# interfaces
export distance, parity_matrix, code_n, code_s, code_k, parity_matrix_x, parity_matrix_z,
rate, metacheck_matrix_x, metacheck_matrix_z, metacheck_matrix, bivariate_bicycle_code_k,
generator_polynomial
export AbstractECC, AbstractQECC, AbstractCECC, AbstractCSSCode, AbstractDistanceAlg

# QEC Codes
export Perfect5, Cleve8, Gottesman

# CSS Codes
export Toric, Bitflip3, Phaseflip3, Shor9, Steane7, Surface, CSS, QuantumReedMuller, Triangular488, Triangular666,
DelfosseReichardt, DelfosseReichardtRepCode, DelfosseReichardt823, QuantumTannerGraphProduct, CyclicQuantumTannerGraphProduct,
TillichZemor, random_TillichZemor_code, BivariateBicycleCodeViaCirculantMat

# Classical Codes
export RepCode, ReedMuller, RecursiveReedMuller, Golay, Hamming, random_Gallager_ldpc, GoppaCode, random_Goppa_code

# utilities
export search_self_orthogonal_rm_codes

include("interface.jl")
include("codes/util.jl")

# Classical Codes
include("codes/classical/hamming.jl")
include("codes/classical/repetition.jl")
include("codes/classical/golay.jl")
include("codes/classical/gallager.jl")
include("codes/classical/goppa.jl")

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
include("codes/quantum/tillichzemor.jl")
include("codes/quantum/generalized_circulant_bivariate_bicycle.jl")

# Reed-Muller Codes
include("codes/classical/reedmuller.jl")
include("codes/classical/recursivereedmuller.jl")
include("codes/quantum/quantumreedmuller.jl")

# Delfosse-Reichardt Codes
include("codes/quantum/delfosse_reichardt_code.jl")
include("codes/quantum/delfosse_reichardt_repcode.jl")
include("codes/quantum/delfosse_reichardt_823_code.jl")

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
