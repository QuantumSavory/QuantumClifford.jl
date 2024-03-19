struct Perfect5 <: AbstractECC end

function iscss(::Type{Perfect5})
    return false
end

parity_checks(c::Perfect5) = S"XZZX_
                               _XZZX
                               X_XZZ
                               ZX_XZ"

distance(c::Perfect5) = 3
