struct SmallestColorCode <: AbstractECC end

code_n(c::SmallestColorCode) = 8

code_k(c::SmallestColorCode) = 8

parity_checks(c::SmallestColorCode) = S"ZZZZIIII
                                        XXXXIIII
                                        IIIIZZZZ
                                        IIIIXXXX
                                        IXYZIXYZ
                                        IZXYIZXY"

distance(c::SmallestColorCode) = 3
