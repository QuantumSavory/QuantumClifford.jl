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


x = parity_checks
#println(x)
y =    S"ZZ_______
         _ZZ______
         ___ZZ____
         ____ZZ___
         ______ZZ_
         _______ZZ
         XXXXXX___
         ___XXXXXX" 
#=                        
#Testing parity check type
typeof(x)
typeof(y)

dim_encondingc = length(y)
println(dim_encondingc)
tracking1 = (dim_encondingc - 1)
println(tracking1)
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
parity_matrix(c::Shor9) = stab_to_gf2(parity_checks(c::Shor9)) #arg included body but not used??
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

#Related functions
#=
code_s(c::Shor9) = (size(parity_checks(c))[1]) / code_n(c)

code_k(c::Shor9) = code_n(c::Shor9) - code_s(c::Shor9)

rate(c::Shor9) = code_k(c::Shor9)/code_s(c::Shor9)
=#
#Naive Syndrome circuit ----------------------------------
naive_sc = []
dim_encondingc = length(y)
println(dim_encondingc)
#= tried:
length(parity_checks((c::Shor9))())
length(parity_checks)
length(parity_checks(Shor9))
length(parity_checks())
=#

ancilla_qubit = dim_encondingc+1
tracking1 = (dim_encondingc - 1)
println(tracking1)
tracking2 = 2

#iterating through all the steps of the encoding circuit
for qubit in 1:tracking1
    #append!(naive_sc, sCNOT(1,ancilla_qubit))
    #append!(naive_sc, sCNOT(tracking2,ancilla_qubit))
    ancilla_qubit + 1
    tracking2 +1
    return naive_sc
end
#=
naive_syndrome_circuit(c::Shor9) = naive_sc
=#
#----------------------------------------------------------------

distance(c::Shor9) = 3 #arg included body but not used

logx_ops(c::Shor9) = P"XXXXXXXXX"
                       
logz_ops(c::Shor9) = P"ZZZZZZZZZ"

isdegenerate(c::Shor9) = true 
