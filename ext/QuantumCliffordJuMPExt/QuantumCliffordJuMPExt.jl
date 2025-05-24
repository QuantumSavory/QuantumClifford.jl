module QuantumCliffordJuMPExt

using DocStringExtensions

import JuMP
import JuMP: @variable, @objective, @constraint, optimize!, Model, set_silent,
    sum, value, optimize!, is_solved_and_feasible, set_time_limit_sec, solution_summary,
    termination_status, MOI.MEMORY_LIMIT, MOI.TIME_LIMIT
import QuantumClifford
import QuantumClifford: stab_to_gf2, logicalxview, logicalzview, canonicalize!,
    MixedDestabilizer, Stabilizer
import QuantumClifford.ECC
import QuantumClifford.ECC: distance, AbstractECC, code_n, code_k, parity_checks,
    logx_ops, logz_ops, parity_checks_x, parity_checks_z, AbstractECC, nqubits,
    DistanceMIPAlgorithm

import SparseArrays
import SparseArrays: SparseMatrixCSC, sparse, spzeros, findnz, sparsevec

import ILog2
import ILog2: ilog2, RoundUp

include("util.jl")
include("min_distance_mixed_integer_programming.jl")

end
