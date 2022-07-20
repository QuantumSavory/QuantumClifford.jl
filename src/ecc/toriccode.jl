#IN PROGRESS
struct Toric <: AbstractECC end 


#Toric x qubit code
parity_checks(c::Toric) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c::Toricx))

syndrome_circuit(c::Toric) = #TODO

#Enconding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#----------------------------------------------------------------

code_n(c::Toric) = #variable input dependent

code_k(c::Toric) = #variable input dependent

code_s(c::Toric) = #variable input dependent

rate(c::Toric) = #variable input dependent

distance(c::Toric) = #variable input dependent

logx_ops(c::Toric) = #variable input dependent

logz_ops(c::Toric) = #variable input dependent

logy_ops(c::Toric) = #variable input dependent

isdegenerate(c::Toric) = #variable input dependent

