using Test
using QuantumClifford
using QuantumClifford.ECC: BCH

@test_throws ArgumentError BCH(2, 2)
@test_throws ArgumentError BCH(3, 4)
