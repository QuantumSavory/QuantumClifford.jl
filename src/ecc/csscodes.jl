import .ECC

#using Nemo 
using LinearAlgebra
using .ECC
#using Statistics

#structure for CSS codes
struct CSS <: AbstractECC 
    H
    G
end

#----------------CSS code generation ----------------------------------------------------------

#-----------Building CSS code ----------------
function parity_checks(c::CSS)
    
    #defining X & Z matrix
    X_matrix = c.H
    Z_matrix = c.G

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

    return Stabilizer(X_bool,Z_bool)

end

code_n(c::CSS) = size(X_bool, 1) #variable input dependant

parity_matrix(c::CSS) = stab_to_gf2(parity_checks(c)) 

#Encoding circuit ----------------------------------

#encoding_circuit(c::CSS) = [] #TODO -> START SYNDROME CIRCUIT
#----------------------------------------------------------------

code_s(c::CSS) = length(parity_checks(c))

code_k(c::CSS) = code_n(c) - code_s(c)

rate(c::CSS) = code_k(c)/code_s(c)

logx_ops(c::CSS) = P"XXXXXXXXX"

logz_ops(c::CSS) = P"ZZZZZZZZZ"

logy_ops(c::CSS) = P"YYYYYYYYY" 