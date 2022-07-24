#IN PROGRESS
struct Toric <: AbstractECC end 


#Toric x qubit code
parity_checks(c::Toric) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c::Toricx))

#Enconding circuit -----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------

#Syndrome circuit -------------------------------------
naive_syndrome_circuit(c::Steane5) = []

#iterating through all the steps of the encoding circuit
for i in encoding_circuit(c::Steane5):
    #iterating through the different physical qubits
    for a in code_n(c::Steane5):
        #second iteration through available physical qubits (for CNOT gates)
        for b in code_n(c::Steane5):
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
code_n(c::Toric) = #variable input dependent

code_k(c::Toric) = #variable input dependent

code_s(c::Toric) = #variable input dependent

rate(c::Toric) = code_k(c::Toric)/code_s(c::Toric)

#distance(c::Toric) = undefined for now

logx_ops(c::Toric) = #variable input dependent

logz_ops(c::Toric) = #variable input dependent

logy_ops(c::Toric) = #variable input dependent

isdegenerate(c::Toric) = #variable input dependent

