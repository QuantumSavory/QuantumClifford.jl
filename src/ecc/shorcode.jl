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

#=                            
#Testing parity check type
x = parity_checks
typeof(x)
=#
##--------------------------------------------------------

#=
parity_checks(c::Shor9) = if c == Shor9
                            S"ZZ_______
                              _ZZ______
                              ___ZZ____
                              ____ZZ___
                              ______ZZ_
                              _______ZZ
                              XXXXXX___
                              ___XXXXXX"
                           end
=#
parity_matrix(c::Shor9) = stab_to_gf2(parity_checks(c::Shor9))
#Enconding circuit ----------------------------------
c1 = sCNOT(1,4)
c2 = sCNOT(1,7)

h1 = sHadamard(1)
h2 = sHadamard(4)
h3 = sHadamard(7)

c3 = sCNOT(1,2)
c4 = sCNOT(4,5)
c5 = sCNOT(7,8)

c6 = sCNOT(1,3)
c7 = sCNOT(4,6)
c8 = sCNOT(7,9) 

encoding_circuit(c::Shor9) = [c1,c2,h1,h2,h3,c3,c4,c5,c6,c7,c8]

#----------------------------------------------------------------

#These functions are current in ECC
#=
naive_syndrome_circuit(c::Shor9) #Syndrome circuit
#steane tutorial, book <- check doc: 

code_s(c::Shor9) = sizeof(parity_checks) / code_n

code_k(c::Shor9) = code_n(c::Shor9) - code_s(c::Shor9)

rate(c::Shor9) = code_k(c::Shor9)/code_s(c::Shor9)
=#

#naive_syndrome_circuit(c::Shor9)
#=
code_s(Shor9)

code_k(Shor9)

rate(Shor9)
=#
distance(c::Shor9) = 3

logx_ops(c::Shor9) = P"XXXXXXXXX"
                       
logz_ops(c::Shor9) = P"ZZZZZZZZZ"

isdegenerate(c::Shor9) = true 
