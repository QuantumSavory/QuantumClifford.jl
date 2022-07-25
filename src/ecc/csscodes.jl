<<<<<<< HEAD
import .ECC

using LinearAlgebra
using .ECC
=======
>>>>>>> 58a4737 (compiling errors)

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

<<<<<<< HEAD
    #ensuring the vector are type Int8
    Xvec = convert(Array{Int8,1}, Xvec)
    Zvec = convert(Array{Int8,1}, Zvec)
=======
code_n(c) = css_n #variable input 
>>>>>>> 58a4737 (compiling errors)


<<<<<<< HEAD
    #resizing the vectors into desired size
    n = size(c.G,1)
    resize!(Xvec, n*n) 
    resize!(Zvec, n*n)
=======
#Encoding circuit ----------------------------------
>>>>>>> 58a4737 (compiling errors)

    #reshappinng X & Z into matrix
    Z = reshape(Zvec, n, n)
    X = reshape(Xvec, n, n)

<<<<<<< HEAD
    #making X & Z into bool
    Z_bool = !=(0).(Z)
    X_bool = !=(0).(X)

    return Stabilizer(X_bool,Z_bool)

end

code_n(c::CSS) = size(X_bool, 1) #variable input dependant
=======
#Syndrome circuit -------------------------------------
function naive_syndrome(encoding_circuit)
    naive_syndrome_circuit = []

    #iterating through all the steps of the encoding circuit
    for i in 1:size(encoding_circuit)
        #iterating through the different physical qubits
        for a in 1:code_n
            #second iteration through available physical qubits (for CNOT gates)
            for b in 1:code_n
                #change qubit order if CNOT gate
                if encoding_circuit[i] == sCNOT(a,b)
                    #adding the steps to the circuit build
                    append!(naive_syndrome_circuit(sCNOT(b,a)))
            
                #Hadamard gates response -> keep step as is
                else
                    append!(naive_syndrome_circuit(encoding_circuit[i]))                
                end
            end
        end
        return naive_syndrome_circuit
    end
end
#----------------------------------------------------------------

code_s(c) = nrow(S)

code_k(c) = css_n - code_s

rate(c) = code_k/code_s

#distance(c::CSS) = undefined for now

logx_ops(c) = P"XXXXXXXXX"

logy_ops(c) = #TODO

logy_ops(c) = P"YYYYYYYYY" 
>>>>>>> 58a4737 (compiling errors)
