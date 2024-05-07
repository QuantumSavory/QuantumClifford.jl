using Test
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: ReedMuller, Goppa

@test_throws ArgumentError ReedMuller(-1, 3)
@test_throws ArgumentError ReedMuller(1, 0)

@test_throws ArgumentError Goppa(2, 1)
@test_throws ArgumentError Goppa(3, 4)
@test_throws ArgumentError Goppa(2, -1)
