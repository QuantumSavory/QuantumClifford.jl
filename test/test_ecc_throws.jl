using Test
using QuantumClifford
using QuantumClifford.ECC: ReedMuller, BCH, ExtendedReedSolomonMDS

@test_throws ArgumentError ReedMuller(-1, 3)
@test_throws ArgumentError ReedMuller(1, 0)

@test_throws ArgumentError BCH(2, 2)
@test_throws ArgumentError BCH(3, 4)

@test_throws ArgumentError ExtendedReedSolomonMDS(2, 1)
@test_throws ArgumentError ExtendedReedSolomonMDS(3, 10)
@test_throws ArgumentError ExtendedReedSolomonMDS(-2, 1)
