#IN PROGRESS
struct Toric <: AbstractECC end 


#Toric x qubit code
parity_checks(c) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

parity_matrix(c) = stab_to_gf2(parity_checks(c))

#Enconding circuit -----------------------------------

encoding_circuit(c) = [] #TODO
#-----------------------------------------------------

#Syndrome circuit -------------------------------------
naive_syndrome_circuit(c) = []

#iterating through all the steps of the encoding circuit
for i in size(encoding_circuit(c)):
    #iterating through the different physical qubits
    for a in code_n(c):
        #second iteration through available physical qubits (for CNOT gates)
        for b in code_n(c):
            #change qubit order if CNOT gate
            if i == sCNOT(a,b):
                #naming the steps
                @eval 
                $(Symbol(:x, i)) = step[$i]
                #adding the steps to the circuit build
                append!(naive_syndrome_circuit(step[$i], sCNOT(b,a)))

            #change X->Z & vice-versage 
            elseif i == Z(a):
                @eval 
                $(Symbol(:x, i)) = step[$i]
                append!(naive_syndrome_circuit(step[$i], X(a)))

            elseif i == X(a):
                @eval 
                $(Symbol(:x, i)) = step[$i]
                append!(naive_syndrome_circuit(step[$i], Z(a)))

            #Hadamard gates response -> keep step as is
            else:
                @eval 
                $(Symbol(:x, i)) = step[$i]
                append!(encoding_circuit(step[$i], i)) 

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

