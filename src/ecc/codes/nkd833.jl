struct NKD833 <: AbstractECC end

parity_checks(c::NKD833) =S"XXXXXXXX
                            ZZZZZZZZ
                            _X_XYZYZ
                            _XZY_XZY
                            _YXZXZ_Y"

