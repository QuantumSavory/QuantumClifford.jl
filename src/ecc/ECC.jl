module ECC

using QuantumClifford

abstract type AbstractECC end

"""The syndrome measurement circuit of a given code."""
function syndrome_circuit end

"""The encoding circuit of a given code."""
function encoding_circuit end

"""The number of physical qubits in a code."""
function code_n end

"""The number of logical qubits in a code."""
function code_k end

"""The number of stabilizer checks in a code."""
function code_s end

"""The rate of a code."""
function rate end

"""The distance of a code."""
function distance end

"""Parity check tableau of a code."""
function parity_checks end

"""Logical X operations of a code."""
function logx_ops end # can be computed from the parity checks

"""Logical Z operations of a code."""
function logz_ops end # can be computed from the parity checks

"""Is the code degenerate"""
function isdegenerate end

# Add others

include("./small_named_codes.jl")
#include("./surface_codes.jl")
#include...
# Shor code
# Steane, 7 qubit, 5 qubit
# repetition codes
# CSS codes
# concatenated codes
# surface and toric codes
# color codes
# LDPC expander and hypergraph codes
# reed muller and reed solomon codes


end
