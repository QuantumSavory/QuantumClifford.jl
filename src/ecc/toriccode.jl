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

#naive_syndrome(encoding_circuit) #Syndrome circuit

code_n(c) = #variable input dependent

code_k(c) = #variable input dependent

code_s(c) = #variable input dependent

rate(c) = code_k(c)/code_s(c)

#distance(c) = undefined for now

logx_ops(c) = #variable input dependent

logz_ops(c) = #variable input dependent

logy_ops(c) = #variable input dependent

isdegenerate(c) = #variable input dependent

