struct N8_K3_D3 <: AbstractECC end

code_n(c::N8_K3_D3) = 8

parity_checks(c::N8_K3_D3) =   S"XXXXXXXX
                                 ZZZZZZZZ
                                 _X_XYZYZ
                                 _XZY_XZY
                                 _YXZXZ_Y"
distance(c::N8_K3_D3) = 3
