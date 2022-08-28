import .ECC
using .ECC

struct Bitflip3 <: AbstractECC end

"""The number of physical qubits in a code."""
code_n(c::Bitflip3) = 3

"""Parity check tableau of a code."""

parity_checks(c::Bitflip3) = S"___
                               _ZZ
                               Z_Z
                               ZZ_"

#Enconding circuit ----------------------------------
c1 = sCNOT(0,1)
c2 = sCNOT(0,3)

encoding_circuit(c::Bitflip3) = [c1,c2]

#----------------------------------------------------------------

logx_ops(c::Bitflip3) = P"XXXXXXXXX"
                       
logz_ops(c::Bitflip3) = P"ZZZZZZZZZ"
