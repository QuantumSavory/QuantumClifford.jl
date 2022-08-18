include("./hammingcodegenerator.jl")

#using Nemo 
using LinearAlgebra
import .ECC
#import .hammingcodegenerator
using .ECC
#using .hammingcodegenerator: H,G

#using Statistics

#structure for CSS codes
struct CSS <: AbstractECC end

#----------------CSS code generation ----------------------------------------------------------
#----------Dual code -------------------
<<<<<<< HEAD
#defining X & Z matrix
X_matrix = H
=======
#Size -not working atm: needs to be fixed
size_row_H = size(H, 1)
size_column_H = size(H, 2)

size_row_G = size(G, 1)
size_column_G = size(G, 2)

#Dual code build
X_zeros = zeros(Int8, size_row_H, size_column_H)
Z_zeros = zeros(Int8, size_row_G, size_column_G)

#Final X & Z matrix
X_matrix = X_zeros
Z_matrix = G

#transforming the matrices into vec
Xvec = vec(X_matrix)
Zvec = vec(Z_matrix)

#ensuring the vector are type Int8
Xvec = convert(Array{Int8,1}, Xvec)
Zvec = convert(Array{Int8,1}, Zvec)


#resizing the vectors into desired size
resize!(Xvec, 7*7)
resize!(Zvec, 7*7)

#reshappinng X & Z into matrix
Z = reshape(Zvec, 7, 7)
X = reshape(Xvec, 7, 7)

#making X & Z into bool
Z_bool = !=(0).(Z)
X_bool = !=(0).(X)

#-----------Building CSS code ----------------
#RE-SIZE BEFORE NOT AFTER
parity_checks(c::CSS) = Stabilizer(X_bool,Z_bool) #READ MANUAL 

code_n(c::CSS) = size(X_bool, 1) #variable input dependant

parity_matrix(c::CSS) = stab_to_gf2(parity_checks) 

#Encoding circuit ----------------------------------

#encoding_circuit(c::CSS) = [] #TODO -> START SYNDROME CIRCUIT
#----------------------------------------------------------------

code_s(c::CSS) = length(parity_checks(c))

code_k(c::CSS) = code_n(c) - code_s(c)

rate(c::CSS) = code_k(c)/code_s(c)


logx_ops(c::CSS) = P"XXXXXXXXX"

logz_ops(c::CSS) = P"ZZZZZZZZZ"

logy_ops(c::CSS) = P"YYYYYYYYY" 

naive_syndrome(c::CSS) #Syndrome circuit
