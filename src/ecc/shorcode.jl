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

parity_matrix(c::Shor9) = stab_to_gf2(parity_checks(c::Shor9))

#Syndrome circuit -------------------------------------
c1 = sCNOT(4,1)
c2 = sCNOT(7,1)

h1 = sHadamard(1) #H change?
h2 = sHadamard(4) #?
h3 = sHadamard(7) #?

c3 = sCNOT(5,4)
c4 = sCNOT(6,4)
c5 = sCNOT(8,7)
c6 = sCNOT(9,7) 

naive_syndrome_circuit(c::Shor9) = [c1,c2,h1,h2,h3,c3,c4,c5,c6]
#----------------------------------------------------------------

#Enconding circuit ----------------------------------
c1 = sCNOT(1,4)
c2 = sCNOT(1,7)

h1 = sHadamard(1)
h2 = sHadamard(4)
h3 = sHadamard(7)

c3 = sCNOT(4,5)
c4 = sCNOT(4,6)
c5 = sCNOT(7,8)
c6 = sCNOT(7,9) 

encoding_circuit(c::Shor9) = [c1,c2,h1,h2,h3,c3,c4,c5,c6]
#----------------------------------------------------------------

code_s(c::Shor9) = nrow(S)

code_k(c::Shor9) = code_n(c::Shor9) - code_s(c::Shor9)

rate(c::Shor9) = code_k(c::Shor9)/code_s(c::Shor9)

distance(c::Shor9) = code_s(c::Shor9)/2 #mmm.... not correct

logx_ops(c::Shor9) = P"XXXXXXXXX"
                       
logz_ops(c::Shor9) = P"ZZZZZZZZZ"

isdegenerate(c::Shor9) = true 
