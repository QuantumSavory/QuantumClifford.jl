"""http://www.codetables.de/QECC.php?q=4&n=5&k=1"""
struct FiveOneThree <: AbstractECC end

code_n(c::FiveOneThree) = 8

parity_checks(c::FiveOneThree) = S"YIZXY
                                    ZXIZY
                                    ZIXYZ
                                    IZZZZ"