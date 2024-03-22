using Test
using QuantumClifford
using QuantumClifford.ECC

@test_throws ArgumentError ReedMuller(-1, 3)
@test_throws ArgumentError ReedMuller(1, 0) 