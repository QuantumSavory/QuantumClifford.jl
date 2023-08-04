"""A pedagogical example of a quantum error correcting [8,3] code used in [cleve1997efficient](@cite)."""
struct Cleve8 <: AbstractECC end

code_n(c::Cleve8) = 8

parity_checks(c::Cleve8) = S"XXXXXXXX
                             ZZZZZZZZ
                             XIXIZYZY
                             XIYZXIYZ
                             XZIYIYXZ"
