import .ECC

using LinearAlgebra
using .ECC

struct Toric <: AbstractECC 
    n #qubits
end 

#Toric n qubit code => nx(n-1) parity check matrix
parity_checks(c::Toric) = S"XZZX_
                            _XZZX
                            X_XZZ
                            ZX_XZ"

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------

code_s(c::Toric) = length(parity_checks(c))

code_k(c::Toric) = c.n - code_s(c)

rate(c::Toric) = code_k(c::Toric)/code_s(c::Toric)

logx_ops(c::Toric) = #variable input dependent

logz_ops(c::Toric) = #variable input dependent

logy_ops(c::Toric) = #variable input dependent


