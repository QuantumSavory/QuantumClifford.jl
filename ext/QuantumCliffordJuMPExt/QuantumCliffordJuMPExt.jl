module QuantumCliffordJuMPExt

using DocStringExtensions

import JuMP
import JuMP: @variable, @objective, @constraint, optimize!, Model, set_silent,
    sum, value
import GLPK
import GLPK: Optimizer

import QuantumClifford
import QuantumClifford: stab_to_gf2, logicalxview, logicalzview, canonicalize!,
    MixedDestabilizer, Stabilizer
import QuantumClifford.ECC
import QuantumClifford.ECC: minimum_distance, AbstractECC, code_n, code_k,
    parity_checks

import SparseArrays
import SparseArrays: SparseMatrixCSC, sparse

include("util.jl")
include("min_distance_mixed_integer_programming.jl")

end
