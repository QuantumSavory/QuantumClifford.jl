using QuantumClifford

@testset "Enumeration" begin
    @test Set([copy(op) for op in enumerate_cliffords(1)]) == Set((random_clifford(1,phases=false) for i in 1:70))
    @test Set([copy(op) for op in enumerate_phases(enumerate_cliffords(1))]) == Set((random_clifford(1,phases=true) for i in 1:300))
    @test length(collect(enumerate_cliffords(2))) == length(collect(enumerate_phases(enumerate_cliffords(2))))/2^4 == 720
    @test first(enumerate_cliffords(3)) == C"X__ __Z _Z_ Z__ __X _X_"
    @test first(enumerate_cliffords(5)) == C"X____ ____Z ___Z_ __Z__ _Z___ Z____ ____X ___X_ __X__ _X___"
end
