# TODO Steane5 qubit code

struct Steane7 <: AbstractECC end

iscss(c::Steane7) = true

parity_checks(c::Steane7) = S"___XXXX
                              _XX__XX
                              X_X_X_X
                              ___ZZZZ
                              _ZZ__ZZ
                              Z_Z_Z_Z"

parity_checks_x(c::Steane7) = stab_to_gf2(parity_checks(Steane7()))[1:3,1:end÷2]
parity_checks_z(c::Steane7) = stab_to_gf2(parity_checks(Steane7()))[4:end,end÷2+1:end]

distance(c::Steane7) = 3
