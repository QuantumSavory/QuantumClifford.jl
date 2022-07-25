module ECC

using QuantumClifford

abstract type AbstractECC end

"""The syndrome measurement circuit of a given code."""
function syndrome_circuits(encoding_circuit) #TODO: add the other 3 types of syndrome circuit
    
    naive_syndrome_circuit = []
    #iterating through all the steps of the encoding circuit
    for i in 1:size(encoding_circuit)
        #iterating through the different physical qubits
        for a in 1:code_n
            #second iteration through available physical qubits (for CNOT gates)
            for b in 1:code_n
                #change qubit order if CNOT gate
                if encoding_circuit[i] == sCNOT(a,b)
                    #adding the steps to the circuit build
                    append!(naive_syndrome_circuit(sCNOT(b,a)))
            
                #Hadamard gates response -> keep step as is
                else
                    append!(naive_syndrome_circuit(encoding_circuit[i]))                        
                end
            end
        end
        return naive_syndrome_circuit
    end

    #fault tolerant (3 types) -Neil, Steane, Shor

end 

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

"""Logical Y operations of a code."""
function logy_ops end # can be computed from the parity checks

"""Is the code degenerate"""
function isdegenerate end

#------IN PROGRESS-------------------------

# Shor code
include("./shorcode.jl")

# Steane, 7 qubit, 5 qubit
include("./steanecode.jl")

# CSS codes
#include("./csscodes.jl")

#------TODO--------------------------------

# Toric codes
#include("./toriccode.jl") #has to b commented for tsting

# Surface codes
#include("./surface_codes.jl")
#include...
# repetition codes
# concatenated codes
# color codes
# LDPC expander and hypergraph codes
# reed muller and reed solomon codes

# Add others ------------------------

end
