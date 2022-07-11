struct Shor9 <: AbstractECC end

parity_checks(c::Shor9) = S"ZZ_ZZ_ZZ_
                            ____ZZ_ZZ
                            Z_ZZZ____
                            ____ZZ___
                            YYXXYY___
                            Z_ZYXYYXY
                             _ZZ_____
                            ____ZZ___"

code_n(c::Shor9) = 9

parity_matrix(c::Shor9) = [[1:0:0:0:0:0:0:0]
                           [1:1:0:0:0:0:0:0]
                           [0:1:0:0:0:0:0:0]
                           [0:0:1:0:0:0:0:0]
                           [0:0:1:1:0:0:0:0]
                           [0:0:0:1:0:0:0:0]
                           [0:0:0:0:1:0:0:0]
                           [0:0:0:0:1:1:0:0]
                           [0:0:0:0:0:1:0:0]]