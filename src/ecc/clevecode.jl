struct Cleve8 <: AbstractECC end

code_n(c::Cleve8) = 8

parity_checks(c::Cleve8) = S"XXXXXXXX
                            ZZZZZZZZ
                            XIXIZYZY
                            XIYZXIYZ
                            XZIYIYXZ"