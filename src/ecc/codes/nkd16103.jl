struct NKD16103 <: AbstractECC end

parity_checks(c::NKD16103) =S"XXXXXXXXXXXXXXXX
                               ZZZZZZZZZZZZZZZZ
                               _X_X_X_XZYZYZYZY 
                               _X_XZYZYX_X_YZYZ
                               _XZYX_YZ_XZYX_YZ
                               _YXZ_YXZ_YXZ_YXZ"