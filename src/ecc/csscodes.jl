import .ECC

using LinearAlgebra
using .ECC

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
    n = size(c.G,1)
    resize!(Xvec, n*n) 
    resize!(Zvec, n*n)

    #reshappinng X & Z into matrix
    Z = reshape(Zvec, n, n)
    X = reshape(Xvec, n, n)

    #making X & Z into bool
    Z_bool = !=(0).(Z)
    X_bool = !=(0).(X)

    return Stabilizer(X_bool,Z_bool)

end

code_n(c::CSS) = size(X_bool, 1) #variable input dependant

#Encoding circuit ----------------------------------

encoding_circuit(c::CSS) = [] #TODO
#----------------------------------------------------------------

logx_ops(c::CSS) = P"XXXXXXXXX"

logz_ops(c::CSS) = P"ZZZZZZZZZ"

logy_ops(c::CSS) = P"YYYYYYYYY" 