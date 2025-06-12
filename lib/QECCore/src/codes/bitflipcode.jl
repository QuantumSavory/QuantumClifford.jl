struct Bitflip3 <: AbstractECC end

code_n(c::Bitflip3) = 3

function iscss(::Type{Bitflip3})
    return true
end

parity_checks(c::Bitflip3) = S"_ZZ
                               Z_Z"
