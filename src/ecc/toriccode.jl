#IN PROGRESS
struct Toric <: AbstractECC end 

#Toric x qubit code
parity_checks(c::Toric) = S"XZZX_
                            _XZZX
                            X_XZZ
                            ZX_XZ"

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------
#These functions are current in ECC
#=
naive_syndrome(c::Toric) #Syndrome circuit

code_n(c::Toric) = #variable input dependent

code_k(c::Toric) = #variable input dependent

code_s(c::Toric) = #variable input dependent

rate(c::Toric) = code_k(c::Toric)/code_s(c::Toric)
=#

logx_ops(c::Toric) = #variable input dependent

logz_ops(c::Toric) = #variable input dependent

logy_ops(c::Toric) = #variable input dependent

isdegenerate(c::Toric) = #variable input dependent

