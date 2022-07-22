struct Steane5 <: AbstractECC end
struct Steane7 <: AbstractECC end


#Steane 5 qubit code
parity_checks(c::Steane5) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

code_n(c::Steane5) = 5

parity_matrix(c::Steane5) = stab_to_gf2(parity_checks(c::Steane5))

#Syndrome circuit -------------------------------------
# Figure 5: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiK5u3rlI35AhVAh_0HHSwhChIQFnoECAgQAQ&url=https%3A%2F%2Fwww.sciencedirect.com%2Ftopics%2Fmathematics%2Fhadamard-gate&usg=AOvVaw0o2EggDP0p5Jni0POFB7u1

c1 = sCNOT(6,1)
c2 = sCNOT(6,2)
C3 = sCNOT(6,3)
c4 = sCNOT(6,4)
C5 = sCNOT(6,5)

x1 = X(4)
x2 = X(5)
x3 = X(6)

naive_syndrome_circuit(c::Steane5) = [c1,c2,c3,c4,c5,x1,x2,x3] 
#----------------------------------------------------------------


#Enconding circuit ----------------------------------
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

code_s(c::Steane5) = nrow(S)

code_k(c::Steane5) = code_n(c::Steane5) - code_s(c::Steane5)

rate(c::Steane5) = code_k(c::Steane5)/code_s(c::Steane5)

distance(c::Steane5) = code_s(c::Steane5)/2

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

#Syndrome circuit -------------------------------------
c1 = sCNOT(2,1)
c2 = sCNOT(3,1)

h1 = sHadamard(5) #do H need to change?
h2 = sHadamard(6) #?
h3 = sHadamard(7) #?

c3 = sCNOT(1,7)
C4 = sCNOT(1,7)
c5 = sCNOT(4,7)
c6 = sCNOT(1,6)
c7 = sCNOT(3,6)
c8 = sCNOT(4,6)
c9 = sCNOT(2,5)
c10 = sCNOT(3,5)
c11 = sCNOT(2,5)

naive_syndrome_circuit(c::Steane7) =  [c1,c2,h1,h2,h3,c3,c4,c5,c6, c7, c8, c9, c10, c11]
#----------------------------------------------------------------

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

encoding_circuit(c::Steane7) = [c1,c2,h1,h2,h3,c3,c4,c5,c6, c7, c8, c9, c10, c11]
#----------------------------------------------------------------

code_s(c::Steane7) = nrow(S)

code_k(c::Steane7) = code_n(c::Steane7) - code_s(c::Steane7)

rate(c::Steane7) = code_k(c::Steane7)/code_s(c::Steane7)

distance(c::Steane7) = code_s(c::Steane7)/2 

logx_ops(c::Steane7) = P"XXXXXXX" 
                       
logz_ops(c::Steane7) = P"ZZZZZZZ" 

isdegenerate(c::Steane7) = false 
