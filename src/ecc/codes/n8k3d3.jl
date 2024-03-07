struct N8K3D3 <: AbstractECC end

parity_checks(c::N8K3D3) = S"XXXXXXXX
                            ZZZZZZZZ
                            _X_XYZYZ
                            _XZY_XZY
                            _YXZXZ_Y"

