# TODO Steane5 qubit code

struct Steane7 <: AbstractECC end

parity_checks(c::Steane7) = S"___XXXX
                              _XX__XX
                              X_X_X_X
                              ___ZZZZ
                              _ZZ__ZZ
                              Z_Z_Z_Z"

code_n(c::Steane7) = 7

distance(c::Steane7) = 3
