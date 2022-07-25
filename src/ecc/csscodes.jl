using Nemo 

#structure for CSS codes
struct CSS <: AbstractECC

    function css_n end
    function classical_code_G_matrix end
    function classical_code_H_matrix end
    function classic_parity_checks end
    function dual_code_prt1 end
    function dual_code_prt2 end

end

#----------------Generating Hamming codes ----------------------------------------------------
#= TODO:
    Number the bits starting from 1: bit 1, 2, 3, 4, 5, 6, 7, etc.
    Write the bit numbers in binary: 1, 10, 11, 100, 101, 110, 111, etc.
    All bit positions that are powers of two (have a single 1 bit in the binary form of their position) are parity bits: 1, 2, 4, 8, etc. (1, 10, 100, 1000)
    All other bit positions, with two or more 1 bits in the binary form of their position, are data bits.
    Each data bit is included in a unique set of 2 or more parity bits, as determined by the binary form of its bit position.
        Parity bit 1 covers all bit positions which have the least significant bit set: bit 1 (the parity bit itself), 3, 5, 7, 9, etc.
        Parity bit 2 covers all bit positions which have the second least significant bit set: bits 2-3, 6-7, 10-11, etc.
        Parity bit 4 covers all bit positions which have the third least significant bit set: bits 4–7, 12–15, 20–23, etc.
        Parity bit 8 covers all bit positions which have the fourth least significant bit set: bits 8–15, 24–31, 40–47, etc.
        In general each parity bit covers all bits where the bitwise AND of the parity position and the bit position is non-zero.
=#

#Testing code for CSS code contruction
#[7,4] Hamming code -> Steane's 7

#classical_code_H_matrix
classical_code_H_matrix(c) = [0:0:1; 0:1:0; 0:1:1; 1:0:0; 1:0:1; 1:1:0; 1:1:1]
#H = classical_code_H_matrix(c)
#ishadamard(H)

#classical_code_G_matrix - dual code
#classical_code_G_matrix(c) = [1:0:0:0; 0:1:0:0; 0:0:1:0; 0:0:0:1] #tst
classical_code_G_matrix(c) = gf2_H_to_G(classical_code_H_matrix)
#----------------------------------------------------------------------------------------------


#----------------CSS code generation ----------------------------------------------------------

#-----------Check matrix ----------------

#create matrix, over space range x
check_matrix_X(c) = classical_code_H_matrix(c)
check_matrix_Z(c) = classical_code_G_matrix(c)

#------------------------------

#-----------Building CSS code ----------------
#=
S = Stabilizer([0x2, 7x3],
                  check_matrix_X,
                  check_matrix_Z)
=#
#=
julia> Stabilizer([0x2, 0x0],
                  Bool[1 1;
                       0 0],
                  Bool[0 0;
                       1 1])
- XX
+ ZZ
=#

parity_checks(c) = S 
print(S) = S #testing

code_n(c) = size(classical_code_H_matrix, 1) + size(classical_code_G_matrix, 1) #variable input 

parity_matrix(c) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c) = []#TODO
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
end
#----------------------------------------------------------------

code_s(c) = nrow(S)

code_k(c) = css_n - code_s

rate(c) = code_k/code_s

#distance(c::CSS) = undefined for now

logx_ops(c) = P"XXXXXXXXX"

logz_ops(c) = P"ZZZZZZZZZ"

logy_ops(c) = P"YYYYYYYYY" 
