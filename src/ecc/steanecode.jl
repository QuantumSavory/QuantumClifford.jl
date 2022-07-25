struct Steane5 <: AbstractECC end
struct Steane7 <: AbstractECC end

#Steane 5 qubit code
parity_checks(c::Steane5) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

code_n(c::Steane5) = 5

#Encoding circuit ----------------------------------
#https://www.researchgate.net/figure/Encoding-circuit-for-the-five-qubit-code-a-Circuit-to-encode-the-logical-minus-state_fig1_337273308

<<<<<<< HEAD
function encoding_circuit(c::Steane5)
    z1 = sZ(1)
    h1 = sHadamard(3)
    h2 = sHadamard(4)
    is1 = sInvPhase(1)
    c1 = sCNOT(3,5)
    c2 = sCNOT(4,2)
    h3 = sHadamard(2)
    c3 = sCNOT(4,5)
    c4 = sCNOT(2,1)
    is2 = sInvPhase(3)
    s2 = sPhase(4)
    is3 = sInvPhase(5)
    s2 = sPhase(1)
    s3 = sPhase(2)
    z2 = sZ(3)
    c5 = sCNOT(1,5)
    h4 = sHadamard(5)
    c6 = sCNOT(2,5)

    return  [z1,h1,h2,is1,c1,c2,h3,c3,c4,is2,s2,is3,s2,s3,z2,c5,h4,c6] 
=======
#Encoding circuit ----------------------------------
c1 = sCNOT(1,6)
c2 = sCNOT(2,6)
C3 = sCNOT(3,6)
c4 = sCNOT(4,6)
C5 = sCNOT(5,6)

z1 = Z(4)
z2 = Z(5)
z3 = Z(5)

encoding_circuit(c::Steane5) = [c1,c2,c3,c4,c5,z1,z2,z3] 
#----------------------------------------------------------------

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
>>>>>>> 58a4737 (compiling errors)
end
#----------------------------------------------------------------

distance(c::Steane5) = 3

isdegenerate(c::Steane5) = false 

##----------------------------------------------------------------------------------------------------------------------------

#steane 7 qubit code
parity_checks(c::Steane7) = S"___XXXX
                              _XX__XX
                              X_X_X_X
                              ___ZZZZ
                              _ZZ__ZZ
                              Z_Z_Z_Z"

code_n(c::Steane7) = 7

#Encoding circuit ----------------------------------
function encoding_circuit(c::Steane7)
    c1 = sCNOT(1,2)
    c2 = sCNOT(1,3)

<<<<<<< HEAD
    h1 = sHadamard(5)
    h2 = sHadamard(6)
    h3 = sHadamard(7)

    c3 = sCNOT(7,1)
    C4 = sCNOT(7,1)
    c5 = sCNOT(7,4)
    c6 = sCNOT(6,1)
    c7 = sCNOT(6,3)
    c8 = sCNOT(6,4)
    c9 = sCNOT(5,2)
    c10 = sCNOT(5,3)
    c11 = sCNOT(5,2)

    return [c1,c2,h1,h2,h3,c3,c4,c5,c6, c7, c8, c9, c10, c11]
=======
#Encoding circuit ----------------------------------
h1 = sHadamard(5)
h2 = sHadamard(6)
h3 = sHadamard(7)

c3 = sCNOT(7,1)
C4 = sCNOT(7,1)
c5 = sCNOT(7,4)
c6 = sCNOT(6,1)
c7 = sCNOT(6,3)
c8 = sCNOT(6,4)
c9 = sCNOT(5,2)
c10 = sCNOT(5,3)
c11 = sCNOT(5,2)

encoding_circuit(c::Steane7) = [c1,c2,h1,h2,h3,c3,c4,c5,c6, c7, c8, c9, c10, c11]
#----------------------------------------------------------------

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
                    append!(naive_syndrome_circuit(encoding_circuit(i)))        
                end
            end
        end
        return naive_syndrome_circuit
    end
>>>>>>> 58a4737 (compiling errors)
end
#----------------------------------------------------------------

distance(c::Steane7) = 3

isdegenerate(c::Steane7) = false 
