module QuantumCliffordJuMPExt

using DocStringExtensions

import JuMP
import JuMP: @variable, @objective, @constraint, optimize!, Model, set_silent,
    sum, value, optimize!, is_solved_and_feasible
import GLPK
import GLPK: Optimizer

import QuantumClifford
import QuantumClifford: stab_to_gf2, logicalxview, logicalzview, canonicalize!,
    MixedDestabilizer, Stabilizer
import QuantumClifford.ECC
import QuantumClifford.ECC: distance, AbstractECC, code_n, code_k, parity_checks

import SparseArrays
import SparseArrays: SparseMatrixCSC, sparse, spzeros, findnz, sparsevec

include("util.jl")
include("min_distance_mixed_integer_programming.jl")

end
