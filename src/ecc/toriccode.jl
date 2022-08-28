import .ECC

using LinearAlgebra
using .ECC

struct Toric <: AbstractECC 
    n #qubits
end 

#Toric n qubit code => nx(n-1) parity check matrix
parity_checks(c::Toric)
    code_n = c.n
    #matrix with 1 zero per conection point
    empty_grid_matrix = zeros(code_n, code_n -1)

    #HOW TO I CONNECT THEM

end

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------
code_n(c::Toric) = c.n


