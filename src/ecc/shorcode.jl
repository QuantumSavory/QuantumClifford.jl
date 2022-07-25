import .ECC
using .ECC

struct Shor9 <: AbstractECC end

"""The number of physical qubits in a code."""
code_n(c::Shor9) = 9

"""Parity check tableau of a code."""

parity_checks(c::Shor9) = S"ZZ_______
                            _ZZ______
                            ___ZZ____
                            ____ZZ___
                            ______ZZ_
                            _______ZZ
                            XXXXXX___
                            ___XXXXXX"

#Enconding circuit ----------------------------------

function encoding_circuit(c::Shor9) 
    c1 = sCNOT(1,4)
    c2 = sCNOT(1,7)

    h1 = sHadamard(1)
    h2 = sHadamard(4)
    h3 = sHadamard(7)

    c3 = sCNOT(1,2)
    c4 = sCNOT(4,5)
    c5 = sCNOT(7,8)

<<<<<<< HEAD
    c6 = sCNOT(1,3)
    c7 = sCNOT(4,6)
    c8 = sCNOT(7,9) 
    
    return [c1,c2,h1,h2,h3,c3,c4,c5,c6,c7,c8]
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
            
                #change X->Z & vice-versage 
                elseif encoding_circuit[i] == Z(a)
                    append!(naive_syndrome_circuit(X(a)))
            
                elseif encoding_circuit[i] == X(a)
                    append!(naive_syndrome_circuit(Z(a)))
            
                #Hadamard gates response -> keep step as is
                else
                    append!(naive_syndrome_circuit(encoding_circuit[i]))                        
                end
            end
        end
        return naive_syndrome_circuit
    end
>>>>>>> 58a4737 (compiling errors)
end
#----------------------------------------------------------------

distance(c::Shor9) = 3 #arg included body but not used

isdegenerate(c::Shor9) = true 
