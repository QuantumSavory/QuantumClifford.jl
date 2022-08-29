import .ECC

using LinearAlgebra
using .ECC

struct Toric <: AbstractECC 
    lx
    lz
end 

code_n(c::Toric) = c.lx * c.lz

#Parity checks ----------------------------------

#matrix with 1 zero per conection point 
empty_grid_matrix(c::Toric) = zeros(c.lx, c.lz)

#HOW TO I CONNECT THEM

parity_checks(c::Toric) = 1 #temporary
#-----------------------------------------------------

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------

distance(c::Toric) = min(c.lx,c.lz)

