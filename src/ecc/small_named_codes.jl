struct Shor9 <: AbstractECC end

parity_checks(c::Shor9) = S"ZZ_______
                            _ZZ______
                            ___ZZ____
                            ____ZZ___
                            ______ZZ_
                            _______ZZ
                            XXXXXX___
                            ___XXXXXX"

code_n(c::Shor9) = 9

#...
