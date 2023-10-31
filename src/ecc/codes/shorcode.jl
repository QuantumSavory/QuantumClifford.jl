struct Shor9 <: AbstractECC end

code_n(c::Shor9) = 9

parity_checks(c::Shor9) = S"ZZ_______
                            _ZZ______
                            ___ZZ____
                            ____ZZ___
                            ______ZZ_
                            _______ZZ
                            XXXXXX___
                            ___XXXXXX"

parity_checks_x(c::Shor9) = stab_to_gf2(parity_checks(Shor9()))[end-1:end,1:end÷2]
parity_checks_z(c::Shor9) = stab_to_gf2(parity_checks(Shor9()))[1:end-2,end÷2+1:end]

distance(c::Shor9) = 3
