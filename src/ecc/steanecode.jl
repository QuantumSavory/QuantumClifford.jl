struct Steane5 <: AbstractECC end
struct Steane7 <: AbstractECC end


#Steane 5 qubit code
parity_checks(c::Steane5) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

code_n(c::Steane5) = 5

parity_matrix(c::Steane5) = stab_to_gf2(parity_checks(c::Steane5))

syndrome_circuit(c::Steane5) = #TODO

#Enconding circuit ----------------------------------
c1 = sCNOT(1,6)
c2 = sCNOT(2,6)
C3 = sCNOT(3,6)
c4 = sCNOT(4,6)
C5 = sCNOT(5,6)
#https://cs269q.stanford.edu/projects2019/stabilizer_code_report_Y.pdf

encoding_circuit(c::Steane5) = [] #TODO: unsure of source #can there be circuits with X&Z log_ops??
#----------------------------------------------------------------

code_k(c::Steane5) = 1

code_s(c::Steane5) = 4 

rate(c::Steane5) = 1/4

distance(c::Steane5) = 3 

logx_ops(c::Steane5) = P"XXXXX" 
                       
logz_ops(c::Steane5) = P"ZZZZZ" 

isdegenerate(c::Steane5) = false 


#steane 7 qubit code
parity_checks(c::Steane7) = S"___XXXX
                              _XX__XX
                              X_X_X_X
                              ___ZZZZ
                              _ZZ__ZZ
                              Z_Z_Z_Z"

code_n(c::Steane7) = 7

parity_matrix(c::Steane7) = stab_to_gf2(parity_checks(c::Steane7))

#start : see p9/31
#start simple <- create funct for all 4 versions
syndrome_circuit(c::Steane7) = #TODO: automate

#Enconding circuit ----------------------------------
c1 = sCNOT(1,2)
c2 = sCNOT(1,3)

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
#TODO: automate
encoding_circuit(c::Steane7) = [c1,c2,h1,h2,h3,c3,c4,c5,c6, c7, c8, c9, c10, c11]
#----------------------------------------------------------------

code_k(c::Steane7) = 1

code_s(c::Steane7) = 6 

rate(c::Steane7) = 1/6 

distance(c::Steane7) = 3 

logx_ops(c::Steane7) = P"XXXXXXX" 
                       
logz_ops(c::Steane7) = P"ZZZZZZZ" 

isdegenerate(c::Steane7) = false 
