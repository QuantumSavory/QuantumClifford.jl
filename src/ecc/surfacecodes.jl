import .ECC

using LinearAlgebra
using .ECC

struct Surface <: AbstractECC 
    lx
    lz
end 

code_n(c::Surface) = c.lx * c.lz

#Parity checks ----------------------------------

#matrix with 1 zero per conection point 
empty_grid_matrix(c::Surface) = zeros(c.lx, c.lz)

#HOW TO I CONNECT THEM

parity_checks(c::Surface) = 1 #temporary
#-----------------------------------------------------

parity_matrix(c::Surface) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Surface) = [] #TODO
#-----------------------------------------------------

distance(c::Surface) = min(c.lx,c.lz)

