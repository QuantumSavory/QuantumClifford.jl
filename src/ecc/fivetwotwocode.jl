"""http://www.codetables.de/QECC.php?q=4&n=5&k=2"""
struct FiveTwoTwo <: AbstractECC end

code_n(c::FiveTwoTwo) = 8

parity_checks(c::FiveTwoTwo) = S"XXXXI
                                IIIIX
                                ZZZZI"