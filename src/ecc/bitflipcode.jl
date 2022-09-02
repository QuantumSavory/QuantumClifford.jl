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
function encoding_circuit(c::Bitflip3)
    c1 = sCNOT(0,1)
    c2 = sCNOT(0,3)

    return [c1,c2]
end 

<<<<<<< HEAD

=======
>>>>>>> a874140 (logx & z ops automation)
