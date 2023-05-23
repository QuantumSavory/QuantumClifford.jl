struct Bitflip3 <: AbstractECC end

code_n(c::Bitflip3) = 3

parity_checks(c::Bitflip3) = S"_ZZ
                               Z_Z"

function encoding_circuit(c::Bitflip3)
    c1 = sCNOT(1,2)
    c2 = sCNOT(1,3)
    return [c1,c2]
end
