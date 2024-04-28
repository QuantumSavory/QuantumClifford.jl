using Test
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: ReedMuller, Goppa

@test_throws ArgumentError ReedMuller(-1, 3)
@test_throws ArgumentError ReedMuller(1, 0)

@test_throws ArgumentError Goppa(6, 1)
@test_throws ArgumentError Goppa(7, 8)
