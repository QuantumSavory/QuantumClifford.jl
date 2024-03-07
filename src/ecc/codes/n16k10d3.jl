struct N16K10D3 <: AbstractECC end

parity_checks(c::N16K10D3) = S"XXXXXXXXXXXXXXXX
                               ZZZZZZZZZZZZZZZZ
                               _X_X_X_XZYZYZYZY 
                               _X_XZYZYX_X_YZYZ
                               _XZYX_YZ_XZYX_YZ
                               _YXZ_YXZ_YXZ_YXZ"