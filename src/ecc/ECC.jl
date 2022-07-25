module ECC

<<<<<<< HEAD
    using QuantumClifford
    abstract type AbstractECC end
    #export AbstractECC

    function sCNOT_gatechange(gate::sCNOT) 
        sCNOT(gate.q2,gate.q1) 
    end

    """The encoding circuit of a given code."""
    function encoding_circuit end 

    """Parity check tableau of a code."""
    function parity_checks  end

    """The number of physical qubits in a code."""
    function code_n end
    
    """The number of stabilizer checks in a code."""
    function code_s(c::AbstractECC) 
        s = length(parity_checks(c))
        return s
    end
    
    """The number of logical qubits in a code."""
    function code_k(c::AbstractECC)
        k = code_n(c) - code_s(c)
        return k
    end
    
    """The rate of a code."""
    function rate(c::AbstractECC)
        rate = code_k(c)/code_s(c)
        return rate
    end
    
    """Naive syndrome circuit""" #TODO: add the other 3 types of syndrome circuit
    function naive_syndrome_circuit(c::AbstractECC)
        naive_sc = []
        dim_encondingc = code_n(c) -1
    
        ancilla_qubit = dim_encondingc
        tracking1 = dim_encondingc 
    
        #iterating through all the steps of the encoding circuit
        for qubit in 0:tracking1
            tracking2 = qubit+1
            if qubit < dim_encondingc
                push!(naive_sc, sCNOT(qubit,ancilla_qubit)) 
                push!(naive_sc, sCNOT(tracking2,ancilla_qubit)) 
                ancilla_qubit = ancilla_qubit + 1
                tracking2 = tracking2 + 1
            end
        end

        return naive_sc
    end

    """The distance of a code."""
    function distance end
    
    """Parity matrix of a code."""
    function parity_matrix(c::AbstractECC) 
        paritym = stab_to_gf2(parity_checks(c::AbstractECC)) 
        return paritym
    end
    
    """Logical X operations of a code."""
    function logx_ops(c::AbstractECC)
        MixedDest = MixedDestabilizer(parity_checks(c))
        logicalxview(MixedDest)
    end 

    """Logical Z operations of a code."""
    function logz_ops(c::AbstractECC)
        MixedDest = MixedDestabilizer(parity_checks(c))
        logicalzview(MixedDest)
    end 

    """Is the code degenerate"""
    function isdegenerate end

    #------IN PROGRESS-------------------------
    # 3 qubit bit flip code 
    include("./bitflipcode.jl")

    # Shor code
    include("./shorcode.jl")

    # Steane, 7 qubit, 5 qubit
    include("./steanecode.jl")

    # CSS codes
    include("./csscodes.jl")

    #Hamming code generator (for CSS codes)
    include("./hammingcodegenerator.jl")

    # Toric codes
    include("./toriccode.jl") 

    # Surface codes
    include("./surfacecodes.jl")

end #module
=======
using QuantumClifford

abstract type AbstractECC end

"""The syndrome measurement circuit of a given code."""
<<<<<<< HEAD
function syndrome_circuit end
=======
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
>>>>>>> 58a4737 (compiling errors)

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
>>>>>>> 27956e1 (initial ECC setup)
