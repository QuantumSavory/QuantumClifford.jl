using Revise 
using BenchmarkTools
using QuantumClifford
using QuantumClifford.ECC: Shor9, naive_syndrome_circuit, code_n, parity_checks, encoding_circuit, code_s, code_k, rate, distance,logx_ops, logz_ops, isdegenerate

"""Inititalizing Shor code"""
code = Shor9()

"""Parity check tableau of a code."""
parity_checks(code)
println("Parity check:", parity_checks(code))

"""The encoding circuit of a given code."""
encoding_circuit(code)
println("Enconding circuit:",encoding_circuit(code))

"""The number of physical qubits in a code."""
code_n(code)
println("Number of physical qubits in a code:",code_n(code))

"""Parity check tableau of a code."""
parity_checks(code)
println("Parity check tableau:",parity_checks(code))

"""The number of stabilizer checks in a code."""
code_s(code) 
println("Stabilizer checks in code:",code_s(code))

"""The rate of a code."""
rate(code)
println("Rate of a code:",rate(code))

"""The distance of a code."""
distance(code)
println("Distance of a code:",rate(code))

"""Logical X operations of a code."""
logx_ops(code)
println("Logical X operations:",logx_ops(code))

"""Logical Z operations of a code."""
logz_ops(code)
println("Logical Z operations:",logz_ops(code))

"""Logical Y operations of a code."""
logy_ops(code)
#println("Logical Y operations:",logy_ops(code))

"""Is the code degenerate"""
isdegenerate(code)
println("Denegerate code:", isdegenerate(code))

##------------------------------------------------------

using Revise 
using BenchmarkTools
using QuantumClifford
using QuantumClifford.ECC: Steane7, naive_syndrome_circuit, code_n, parity_checks, encoding_circuit, code_s, code_k, rate, distance,logx_ops, logz_ops, isdegenerate

"""Inititalizing Shor code"""
code = Steane7()

"""Parity check tableau of a code."""
parity_checks(code)
println("Parity check: \n", parity_checks(code))

"""The encoding circuit of a given code."""
encoding_circuit(code)
println("Enconding circuit:",encoding_circuit(code))

"""The number of physical qubits in a code."""
code_n(code)
println("Number of physical qubits in a code:",code_n(code))

"""Parity check tableau of a code."""
parity_checks(code)
println("Parity check tableau:",parity_checks(code))

"""The number of stabilizer checks in a code."""
code_s(code) 
println("Stabilizer checks in code:",code_s(code))

"""The rate of a code."""
rate(code)
println("Rate of a code:",rate(code))

"""The distance of a code."""
distance(code)
println("Distance of a code:",rate(code))

"""Logical X operations of a code."""
logx_ops(code)
println("Logical X operations:",logx_ops(code))

"""Logical Z operations of a code."""
logz_ops(code)
println("Logical Z operations:",logz_ops(code))

"""Logical Y operations of a code."""
logy_ops(code)
#println("Logical Y operations:",logy_ops(code))

"""Is the code degenerate"""
isdegenerate(code)
println("Denegerate code:", isdegenerate(code))