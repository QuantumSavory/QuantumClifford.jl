using Test
using QuantumClifford
using QuantumClifford.ECC: SteaneReedMuller, ReedMuller

@test_throws ArgumentError ReedMuller(-1, 3)
@test_throws ArgumentError ReedMuller(1, 0) 


@test_throws ArgumentError SteaneReedMuller(-1, 3)
@test_throws ArgumentError SteaneReedMuller(1, 0) 