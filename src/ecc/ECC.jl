module ECC

using QuantumClifford

abstract type AbstractECC end

function sCNOT_gatechange(gate::sCNOT) sCNOT(gate.q2,gate.q1) end


function naive_syndrome_circuit(encondingc::Vector) #TODO: add the other 3 types of syndrome circuit
    
    naive_sc = []
    #iterating through all the steps of the encoding circuit
    dim_encondingc = length(encondingc)
    #=for h in encondingc
        dim_encondingc = dim_encondingc + 1
    end=#

    for i in 1:dim_encondingc
        #change qubit order if CNOT gate
        if encondingc[i] == sCNOT
            #adding the steps to the circuit build
            append!(naive_sc(sCNOT_gatechange(encondingc[i])))
            
        #Hadamard gates response -> keep step as is
        else
            append!(naive_sc(encondingc[i]))                        
        end
    end

    return naive_sc
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
include("./csscodes.jl")

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
