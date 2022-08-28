import .ECC

using LinearAlgebra
using .ECC

struct Toric <: AbstractECC 
    n #qubits
end 

code_n(c::Toric) = c.n

#Parity checks ----------------------------------

#matrix with 1 zero per conection point 
empty_grid_matrix = zeros(code_n, code_n -1)
    
#HOW TO I CONNECT THEM

parity_checks(c::Toric) = 1 #temporary
#-----------------------------------------------------

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------


