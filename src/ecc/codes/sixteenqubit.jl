struct SixteenQubit <: AbstractECC end

code_n(c::SixteenQubit) = 16

parity_checks(c: SixteenQubit) = S"XXXXXXXXXXXXXXXX
                                   ZZZZZZZZZZZZZZZZ
                                   _X_X_X_XZYZYZYZY 
                                   _X_XZYZYX_X_YZYZ
                                   _XZYX_YZ_XZYX_YZ
                                   _YXZ_YXZ_YXZ_YXZ"
distance(c::SixteenQubit) = 3
