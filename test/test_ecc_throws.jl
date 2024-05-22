using Test
using QuantumClifford
using QuantumClifford.ECC: ReedMuller, BCH

@test_throws ArgumentError ReedMuller(-1, 3)
@test_throws ArgumentError ReedMuller(1, 0)

@test_throws ArgumentError BCH(2, 2, 2)
@test_throws ArgumentError BCH(2, 3, 4)
@test_throws ArgumentError BCH(4, 3, 4)
