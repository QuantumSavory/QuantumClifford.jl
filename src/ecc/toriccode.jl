#IN PROGRESS
struct Toric <: AbstractECC end 


#Toric x qubit code
parity_checks(c) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

parity_matrix(c) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c) = [] #TODO
#-----------------------------------------------------

#Syndrome circuit -------------------------------------
function naive_syndrome(encoding_circuit)
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
end
#----------------------------------------------------------------

code_n(c) = #variable input dependent

code_k(c) = #variable input dependent

code_s(c) = #variable input dependent

rate(c) = code_k(c)/code_s(c)

#distance(c) = undefined for now

logx_ops(c) = #variable input dependent

logz_ops(c) = #variable input dependent

logy_ops(c) = #variable input dependent

isdegenerate(c) = #variable input dependent

