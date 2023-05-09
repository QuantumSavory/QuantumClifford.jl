using QuantumClifford

using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
using Test

@testset "Stabilizer to Clifford Conversion" begin
    for i in 1:10
        stab = one(Stabilizer, 10)        
        cl = CliffordOperator(stab)
        cltab = tab(cl)
        @test Stabilizer(cltab[end÷2+1:end]) == stab
    end
    for i in 1:10
        stab = random_stabilizer(10)
        gnd = one(stab)
        cl = CliffordOperator(stab)
        @test apply!(gnd, cl) == stab
    end
    ## TODO: More tests??
end

@testset "Basic Stabs to Clifford" begin
    stabSWAP = S"IZ
                 ZI"
    @test CliffordOperator(stabSWAP) == tSWAP
    ## TODO: More tests??
end
