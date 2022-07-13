struct Steane5 <: AbstractECC end
struct Steane7 <: AbstractECC end


#Steane 5 qubit code
parity_checks(c::Steane5) = S"XZZX_
                              _XZZX
                              X_XZZ
                              ZX_XZ"

code_n(c::Steane5) = 5

parity_matrix(c::Steane5) = stab_to_gf2(parity_checks(c::Steane5))


#steane 7 qubit code
parity_checks(c::Steane7) = S"___XXXX
                              _XX__XX
                              X_X_X_X
                              ___ZZZZ
                              _ZZ__ZZ
                              Z_Z_Z_Z"

code_n(c::Steane7) = 7

parity_matrix(c::Steane7) = stab_to_gf2(parity_checks(c::Steane7))


