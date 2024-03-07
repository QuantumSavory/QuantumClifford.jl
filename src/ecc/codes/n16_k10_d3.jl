struct N16_K10_D3 <: AbstractECC end

code_n(c::N16_K10_D3) = 16

parity_checks(c::N16_K10_D3) =   S"XXXXXXXXXXXXXXXX
                                   ZZZZZZZZZZZZZZZZ
                                   _X_X_X_XZYZYZYZY 
                                   _X_XZYZYX_X_YZYZ
                                   _XZYX_YZ_XZYX_YZ
                                   _YXZ_YXZ_YXZ_YXZ"
distance(c::N16_K10_D3) = 3
