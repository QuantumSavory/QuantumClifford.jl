struct EightQubit <: AbstractECC end

code_n(c::EightQubit) = 8

parity_checks(c::EightQubit) = S"XXXXXXXX
                               ZZZZZZZZ
                               _X_XYZYZ
                               _XZY_XZY
                               _YXZXZ_Y"
distance(c::EightQubit) = 3
