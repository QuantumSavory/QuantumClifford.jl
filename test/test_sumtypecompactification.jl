using Test
using QuantumClifford

@testset "SumTypes compactification" begin
    @test_warn "Could not compactify the circuit" QuantumClifford.pftrajectories([ClassicalXOR{17}((65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81), 282)])
    QuantumClifford.compactify_circuit([ClassicalXOR{3}((65, 66, 67), 282)])
end
