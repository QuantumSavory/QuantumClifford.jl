# TODO [[4,2,2]] qubit code (Simplest Surface Code)

struct code422 <: AbstractECC end

parity_checks(c::code422) = S"XZZX
                              YXXY"
#not sure about the two lines,I have a hunch it's correct. Need verification.
parity_checks_x(c::code422) = stab_to_gf2(parity_checks(code422()))[1:2,1:end÷2]
parity_checks_z(c::code422) = stab_to_gf2(parity_checks(code422()))[2:end,end÷2+1:end]

distance(c::code422) = 2
