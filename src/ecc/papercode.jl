struct Paper8 <: AbstractECC end

code_n(c::Paper8) = 8

parity_checks(c::Paper8) = S"XXXXXXXX
                            ZZZZZZZZ
                            XIXIZYZY
                            XIYZXIYZ
                            XZIYIYXZ"