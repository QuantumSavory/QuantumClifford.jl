using Test
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: ReedMuller, ReedSolomon

@test_throws ArgumentError ReedMuller(-1, 3)
@test_throws ArgumentError ReedMuller(1, 0)

@test_throws ArgumentError ReedSolomon(31, 1)
@test_throws ArgumentError ReedSolomon(31, 35)
@test_throws ArgumentError ReedSolomon(-1, 1)
